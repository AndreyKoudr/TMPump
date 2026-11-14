/*
BSD 2-Clause License

Copyright (c) 2020, Andrey Kudryavtsev (andrewkoudr@hotmail.com)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "TCAD.h"

#include "Conforming3D.h"
#include "OOctree.h"
#include "Particle.h"
#include "OctreeAndTriangles.h"

#include <ctime>
#include <random>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iostream>
#include <limits>
#include <chrono>
#include <assert.h>

/*******************************************************************************
  High-vacuum flow inside one-stage turbomolecular pump with DSMC.
*******************************************************************************/

/**** INPUT *******************************************************************/

// save debug info
//#define SAVE_DEBUG

// blade, located in the project folder; Z-axis is along span, X - along chord
//#define BLADE_FILENAME "NACA.stl"
#define BLADE_FILENAME "TM pump stator.001.STL"

// octree parameters
#define MAX_LEVEL 10

// raster cells will refined to max REFINE_MAX, very important!
// we use uniform background cells for particle positions, but refine level > 1
// to calculate volume, classify inside/outside cells and find inner place 
// places for original particles
#define REFINE_MAX 2

// blade represented by signed distance field
#define BLADE_OUTFILENAME "NACA_SDF.stl"

/**** INPUT *******************************************************************/

// all about doubles
#define T double

// times to decimate points on every edge; decrease it to make triangle
// mesh more refined
#define DECIMATE_STATOR 20
#define DECIMATE_ROTOR 20

// increment list of particles by this value
#define PARTICLES_INCREMENT 1024

// temperature, K
const T temperature = 293.15; 
const T wall_temperature = 293.15; 

// dry air kg/mol
const T molar_mass = 0.02896;

// particle speed along Z at entry plane
const T entry_speed_z = 0.0;
//const T entry_speed_z = -1.0;

// pump radius
const T R = 0.12;

// shaft radius
const T Rshaft = 0.05;

// distance between blades of rotor and stator to define Kn number
const T d = 0.0007;

// number of rotor blades
int rotor_numblades = 16;
// number of stator blades per one rotor blade (per period)
int stator_blades_per_period = 3;

// rotations per minute
const T RPM = 10000.0;

// target P, Pa (vacuum)
const T Ptarget = 1.0e-3;

// rotor/stator blade angles of attack, degrees, supplied to the function with
// opposite sign; the rotor rotates around Z counter-clockwise
const T rotor_angle_attack = 45.0;
const T stator_angle_attack = -30.0;

// cell size
const T cell_size = R / 100.0;

// # particles per cell
const int particles_per_cell = 10;

// stop iteration time, ~10^7 steps for 1 second
const T t_stop = 0.5;

/*****************************************************************************/

using namespace tcad;
using namespace std;
using namespace std::chrono;

int main()
{
//  // particles
//  std::vector<Particle<T>> particles;
//
//#ifdef SAVE_DEBUG
//  {
//    TTriangles<T> tris;
//    makeParticleTriangles<T>(particles,tris,0.0012 * 0.1);
//    tris.saveSTL(DEBUG_DIR + "particles.stl","Particles",true);
//  }
//#endif

  // random generator
  std::mt19937 gen((unsigned int) time(nullptr));

  // start hight-res timer
  high_resolution_clock::time_point t1,t2;
  t1 = high_resolution_clock::now();

  // rotor and stator as trimmed surfaces
  TBrep<T> rotor,stator;

  cout << "Making TM pump rotor and stator (takes time) ..." << endl << endl;

  bool ok = makeTMPump(R,Rshaft,true, // if false, non-periodic is slow due to multiple boolean operations when adding blades
    d,
    //0.062, // ~= PI20 * R / rotor_numblades
    //0.013, // ~= PI20 * Rshaft / stator_numblades
    -1.0,-1.0, // let it set chords itself
    rotor,stator,
    rotor_numblades,stator_blades_per_period,
    TMPAirfoil<T>,TMPContour<T>,
    TMPAirfoil<T>,TMPContour<T>,
    -rotor_angle_attack,-stator_angle_attack);

  ASSERT(ok);

#if 0
  // period for rotor
  T period = 360.0 / rotor_numblades;

  // make second blade of rotor rotated clockwise around Z by period;
  // rotor1 is trailing after rotor
  TBrep<T> rotor1 = rotor;

  TTransform<T> t;
  t.Rotate(TPoint<T>(0.0,0.0,1.0),-period * CPI);
  rotor1.makeTransform(&t);

  // do not use "+" for simplicity
  rotor.addFaces(rotor1);
#endif

  cout << "Making triangles from trimmed surfaces ..." << endl << endl;

  // rotor and stator triangles
  TTriangles<OREAL> rotortris,statortris;

  // boundary points decimated ... times
  ok = rotor.makeTriangles(rotortris,DECIMATE_ROTOR,-1.0,0.1,true) && 
    stator.makeTriangles(statortris,DECIMATE_STATOR,-1.0,0.1,true);

  ASSERT(ok);

#ifdef SAVE_DEBUG
  {
    rotortris.saveSTL(DEBUG_DIR + "rotor.stl","Rotor",true);
    statortris.saveSTL(DEBUG_DIR + "stator.stl","Stator",true);
  }
#endif

  // rotor/stator surfaces grouped by surface type
  std::map<int,std::vector<int>> rotor_faces;
  std::map<int,std::vector<int>> stator_faces;

  for (int i = 0; i < int(rotor.surfaces.size()); i++)
  {
    rotor_faces[rotor.surfaces[i]->type].push_back(i);
  }

  for (int i = 0; i < int(stator.surfaces.size()); i++)
  {
    stator_faces[stator.surfaces[i]->type].push_back(i);
  }

  // rotor/stator triangles grouped by surface type
  std::map<int,TTriangles<T>> rotor_tris;
  std::map<int,TTriangles<T>> stator_tris;

  // rotor : create faces for every face type
  groupTriangles<T>(rotor,rotor_faces,"rotor",DECIMATE_ROTOR,rotor_tris);

  // stator : group triangles for every face type
  groupTriangles<T>(stator,stator_faces,"stator",DECIMATE_STATOR,stator_tris);

  //===== SOLUTION =====

  // these are important parameters for octree refinements

  // molecular mass
  T molecular_mass = m(molar_mass);

  // time step from cell size
  T dt0 = dtFromCellSize<T>(cell_size,temperature,molecular_mass);
  // time step from d - distance between rotor and stator
  T max_blade_speed = maxBladeSpeed(R,RPM);
  T dt1 = dtFromBladeDistance(d,max_blade_speed);
  // actual time step
  T dt = std::min<T>(dt0,dt1) * 0.5;

  // mean thermal speed
  T vmean = meanThermalSpeed(temperature,molecular_mass);
  T vmax = maxThermalSpeed(temperature,molecular_mass);

  // max blade speed
  T vbmax = maxBladeSpeed<T>(R,RPM);

  // max distance a particle can cover during one time step
  T max_step = std::max<T>(vmax,vbmax) * dt;

  // calculate refine level limited by REFINE_MAX
  int refine_level = int(log2(cell_size / max_step));
  LIMIT(refine_level,1,REFINE_MAX); // min 1 as we always need raster cells

  cout << "Making rotor octree with refinement level " << refine_level << "..." << endl << endl;

  // check input parameters!
  OOctree rotree; 

  int inindexrotor = 1;
  int outindexrotor = 2;
  makeOctreeFromTriangles<T>(rotortris,rotree,cell_size,refine_level,
    inindexrotor,outindexrotor,MAX_LEVEL,1.2);

  // calculate cuboid and create octree
  std::pair<OVECTOR,OVECTOR> minmax = statortris.minmax();
  //!!!!!!! this is approximate, calculate volume by inner cells
  TPoint<T> dom = minmax.second - minmax.first;

  cout << endl << "Making stator octree with refinement level " << refine_level << "..." << endl << endl;

  // check input parameters!
  OOctree otree; 

  int inindexstator = 1;
  int outindexstator = 2;
  makeOctreeFromTriangles<T>(statortris,otree,cell_size,refine_level,
    inindexstator,outindexstator,MAX_LEVEL,1.1);

  auto Isize = otree.IJKnumCells.X;
  auto Jsize = otree.IJKnumCells.Y;
  auto Ksize = otree.IJKnumCells.Z;
  int num_cells = int(Isize * Jsize * Ksize);

  // approximate volume
  T volume = dom.X * dom.Y * dom.Z;
  T inner_volume = otree.getVolume(otree.cells(),inindexstator);

  cout << "Domain : " << Isize << "x"  << Jsize << "x"  << Ksize << " = "  << num_cells << 
    " volume " << volume << " inner volume " << inner_volume << endl << endl;

  // # molecules in m^3
  T n = numberDensity<T>(temperature,Ptarget);
  // total number of molecules inside simulation volume
  T Np = Nphysical<T>(temperature,volume,Ptarget);
  // specific weight Fnum (number of molecules in one particle)
  T fnum = Fnum<T>(temperature,volume,num_cells,Ptarget,particles_per_cell); 

  // only sector
  T inlet_area = PI10 * R * R / T(rotor_numblades);
  int nstep = Nstep<T>(n,inlet_area,dt,vmean,fnum);

  // total number of particles
  int num_particles = num_cells * particles_per_cell;

  cout << "Parameters :" << endl;
  cout << "basic cell size            " << cell_size << endl;
  cout << "molecular mass m           " << molecular_mass << endl;
  cout << "number density n           " << n << endl;
  cout << "total # molecules          " << Np << endl;
  cout << "specific weight Fnum       " << fnum << endl;
  cout << "time step from cell size   " << dt0 << endl;
  cout << "time step from d           " << dt1 << endl;
  cout << "time step                  " << dt << endl;
  cout << "mean thermal speed Vn      " << vmean << endl;
  cout << "max thermal speed          " << vmax << endl;
  cout << "max blade speed            " << vbmax << endl;
  cout << "max distance per dt        " << max_step << endl;
  cout << "Nstep                      " << nstep << endl;
  cout << "total number of particles  " << num_particles << endl;
  cout << endl;

  //===== OCTREE =====

  OREAL cellSize = otree.meanRealCellSize(0); // we use background cells
  OREAL minSizeRotor = otree.minRealCellSize(refine_level);
  OREAL minSizeStator = otree.minRealCellSize(refine_level);

  cout << "Octrees :" << endl;
  cout << "number of cells            " << Isize * Jsize * Ksize << endl;
  cout << "background cell size       " << cellSize << endl;
  cout << "rotor min cell size        " << minSizeRotor << endl;
  cout << "stator min cell size       " << minSizeStator << endl;
  cout << endl;

  // save background cells of octree
  #ifdef SAVE_DEBUG
  {
    TTriangles<OREAL> tris;
    otree.cellsToTriangles(otree.cells(),tris);
    tris.saveSTL(DEBUG_DIR + "originalcells.stl","Octree",true);
  }
  #endif

  cout << "Filling octree cells with triangles..." << endl << endl;

  // fill octree cells with STATIC triangles
  //       cell ID        list of...       FaceType and triangle
  std::map<IPosition,std::vector<std::pair<int,std::array<TPoint<T>,3>>>,IPosCompare> stator_cell_tris;
  distributeTrisOverCells<T>(stator_tris,otree,refine_level,3,minSizeStator * 0.5,stator_cell_tris);

  std::map<IPosition,std::vector<std::pair<int,std::array<TPoint<T>,3>>>,IPosCompare> rotor_cell_tris;
  distributeTrisOverCells<T>(rotor_tris,otree,refine_level,3,minSizeRotor * 0.5,rotor_cell_tris);

  //===== PARTICLES =====

  cout << "Filling inner cells with particles..." << endl << endl;

  // particles
  std::vector<Particle<T>> particles;

  // all these particles are set active
  fillInnerCells<T>(otree,inindexstator,rotree,inindexrotor,refine_level,particles,particles_per_cell,gen);
  // ... e.g. #inner cells 687299 #particles 369647

  // set original velocity for all particles inside domain
  for (auto &p : particles)
  {
    setParticleVelocity<T>(p,entry_speed_z,gen,temperature,molecular_mass);
  }

  // only stator, not very meaningful
  LINT Ninner = otree.countCells(otree.cells(),inindexstator);
  cout << "# inner cells " << Ninner << " # particles " << particles.size() << endl << endl;

#ifdef SAVE_DEBUG
  {
    TTriangles<T> tris;
    makeParticleTriangles<T>(particles,tris,cellSize * 0.01,10);
    tris.saveSTL(DEBUG_DIR + "particles.stl","Particles",true);
  }
#endif

  cout << "Calculation :" << endl << endl;

  // period for rotor
  T period = 360.0 / rotor_numblades;

  // time 
  T t = 0.0;
  int step = 0;
  while(t < t_stop)
  {
 //   cout << "\rtime " << t << " sec, step " << step << " # particles " << countParticles<T>(particles) << std::flush;
    cout << "\rtime " << t << " sec, step " << step << std::flush;

    // inject new Nstep particles at entry plane
    for (int i = 0; i < nstep; i++)
    {
      // inject Nstep new particles at entry (CYLINDER_TOP) plane
      int index = newParticleIndex<T>(particles,PARTICLES_INCREMENT);

      // place new particle on top
      setParticleEntryPosition<T>(particles[index],stator_tris[CYLINDER_TOP],gen);

      // set velocity
      setParticleVelocity<T>(particles[index],entry_speed_z,gen,temperature,molecular_mass);

      // set active
      particles[index].active = true;
    }

    // move particles
    for (int i = 0; i < int(particles.size()); i++) 
    {
      if (!particles[i].active) 
        continue;

      // move particle continuously in 3D space
      TPoint<T> oldpos = particles[i].position;
      TPoint<T> newpos = oldpos + particles[i].velocity * dt;
      particles[i].position = newpos;

      // check stator
      bool found = false;
      IPosition scentre;
      if (otree.findCellCentre(newpos,refine_level,scentre))
      {
        // raster
        if (otree.cells()[scentre].data.first == 3)
        {
#if 0
          // find all neighbours
          std::vector<IPosition> centres;
          otree.getAllSameSizeNeighbours(scentre,centres);

          // test all neighbours
          for (auto &ncentre : centres)
          {
            // std::map<IPosition,std::vector<std::pair<int,std::array<TPoint<T>,3>>>,IPosCompare> stator_cell_tris
            if (stator_cell_tris.count(ncentre))
            {
              for (auto &list : stator_cell_tris[ncentre])
              {
                // type
                auto type = list.first;
                T U = 0.0;
                TPoint<T> intr;
                bool backfaced = false;
                if (segTriIntersectMT<T>(oldpos,newpos,list.second/* tri */,U,intr,PARM_TOLERANCE,&backfaced))
                {
                  if (backfaced)
                  {
                    // this must not be, exclude this particle
                    particles[i].active = false;
                  } else
                  {
                    // out of chamber through exit or entry
                    if (type == CYLINDER_TOP || type == CYLINDER_BOTTOM)
                    {
                      particles[i].active = false;
                    }
                  }
                  found = true;
                  break;
                }
              }
            }
            if (found)
              break;
          }

#else
          // test all tris
          for (auto &tris : stator_tris)
          {
            T U = 0.0;
            TPoint<T> intersection,normal;
            bool backfaced = false;
            int index = tris.second.segIntersectAll(oldpos,newpos,U,intersection,normal,backfaced,PARM_TOLERANCE);
            if (index >= 0)
            {
              if (backfaced)
              {
                // this must not be, exclude this particle
                particles[i].active = false;
              } else
              {
                // out of chamber through exit or entry
                if (tris.first == CYLINDER_TOP || tris.first == CYLINDER_BOTTOM)
                {
                  particles[i].active = false;
                } else
                {
                  // reflect from stationary wall (all stator surfaces)
                  //TPoint<T> wall_velocity;
                  //applyParticleDiffuseReflection<T>(particles[i], 
                  //  wall_velocity,normal,wall_temperature,molecular_mass,gen);
                }
              }
            }
          }
#endif
        }
      }
    }

    // rotate rotor

    t += dt;
    step++;
  }

  cout << endl << "Press ENTER to exit..." << endl;
  getchar();
}
