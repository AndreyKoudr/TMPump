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
#define SAVE_DEBUG

// blade, located in the project folder; Z-axis is along span, X - along chord
//#define BLADE_FILENAME "NACA.stl"
#define BLADE_FILENAME "TM pump stator.001.STL"

// create an octree?
//#define USE_OCTREE

// octree parameters
#define IMAX 4
#define JMAX 3
#define KMAX 2
#define MAX_LEVEL 10

// raster cells will refined to REFINE_MAX, very important!
#define REFINE_MAX 6

// blade represented by signed distance field
#define BLADE_OUTFILENAME "NACA_SDF.stl"

/**** INPUT *******************************************************************/

// all about doubles
#define T double

// times to decimate points on every edge
#define DECIMATE 20

// temperature, K
const T temperature = 293.15; 

// dry air kg/mol
const T molar_mass = 0.02896;

// particle speed along Z at entry plane
const T entry_speed_z = -1.0;

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
const T RPM = 10000;

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

/*****************************************************************************/

using namespace tcad;
using namespace std;
using namespace std::chrono;

int main()
{
  // random generator
  std::mt19937 gen((unsigned int) time(nullptr));

  // start hight-res timer
  high_resolution_clock::time_point t1,t2;
  t1 = high_resolution_clock::now();

  // rotor and stator as trimmed surfaces
  TBrep<T> rotor,stator;

  // rotor and stator triangles
  TTriangles<OREAL> rotor_tris,stator_tris;

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

  cout << "Making triangles from trimmed surfaces ..." << endl << endl;

  // boundary points decimated ... times
  ok = rotor.makeTriangles(rotor_tris,DECIMATE,-1.0,0.1,true) && 
    stator.makeTriangles(stator_tris,DECIMATE,-1.0,0.1,true);

  ASSERT(ok);

  rotor_tris.saveSTL(DEBUG_DIR + "rotor.stl","Rotor",true);
  stator_tris.saveSTL(DEBUG_DIR + "stator.stl","Stator",true);

  // rotor/stator faces grouped by type
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

  // rotor : create faces for every face type
  for (auto t : rotor_faces)
  {
    TTriangles<T> typetris;
    for (int i = 0; i < int(t.second.size()); i++)
    {
      TTriangles<T> tris;
      rotor.makeFaceTriangles(t.second[i],tris,DECIMATE,-1.0,0.1,true);
      typetris.addTriangles(tris);
    }
    typetris.saveSTL(DEBUG_DIR + "rotor_faces_" + facetypeNames[t.first] + ".stl","Rotor",true);
  }

  // stator : create faces for every face type
  TTriangles<T> entry_tris;
  for (auto t : stator_faces)
  {
    TTriangles<T> typetris;
    for (int i = 0; i < int(t.second.size()); i++)
    {
      TTriangles<T> tris;
      stator.makeFaceTriangles(t.second[i],tris,DECIMATE,-1.0,0.1,true);
      typetris.addTriangles(tris);
    }
    typetris.saveSTL(DEBUG_DIR + "stator_faces_" + facetypeNames[t.first] + ".stl","Stator",true);
    if (t.first == CYLINDER_TOP)
    {
      entry_tris = typetris;
    }
  }

  //===== SOLUTION =====

  // calculate cuboid and create octree
  std::pair<OVECTOR,OVECTOR> minmax = stator_tris.minmax();

  OVECTOR min = minmax.first;
  OVECTOR max = minmax.second;
  OVECTOR dim = (max - min) * OREAL(0.5);
  OVECTOR c = (max + min) * OREAL(0.5);
  OVECTOR omin = c - dim * OREAL(1.1);
  OVECTOR omax = c + dim * OREAL(1.1);

  // calculate domain size in cells
  OVECTOR dom = omax - omin;
  int Isize = ROUND(dom.X / cell_size); 
  LIMIT_MIN(Isize,1);
  int Jsize = ROUND(dom.Y / cell_size); 
  LIMIT_MIN(Jsize,1);
  int Ksize = ROUND(dom.Z / cell_size); 
  LIMIT_MIN(Ksize,1);
  int num_cells = Isize * Jsize * Ksize;

  cout << "Domain : " << Isize << "x"  << Jsize << "x"  << Ksize << " = "  << num_cells << endl << endl;

  // molecular mass
  T molecular_mass = m(molar_mass);

  // # molecules in m^3
  T n = numberDensity<T>(temperature,Ptarget);
  
  // approximate volume
  T volume = dom.X * dom.Y * dom.Z;
  // total number of molecules inside simulation volume
  T Np = Nphysical<T>(temperature,volume,Ptarget);
  // specific weight Fnum (number of molecules in one particle)
  T fnum = Fnum<T>(temperature,volume,num_cells,Ptarget,particles_per_cell); 

  // time step from cell size
  T dt0 = dtFromCellSize<T>(cell_size,temperature,molecular_mass);
  // time step from d - distance between rotor and stator
  T max_blade_speed = maxBladeSpeed(R,RPM);
  T dt1 = dtFromBladeDistance(d,max_blade_speed);
  // actual time step
  T dt = std::min<T>(dt0,dt1) * 0.5;

  // mean thermal speed
  T vn = meanThermalSpeed(temperature,molecular_mass);

  // only sector
  T inlet_area = PI10 * R * R / T(rotor_numblades);
  int nstep = Nstep<T>(n,inlet_area,dt,vn,fnum);

  // total number of particles
  int num_particles = num_cells * particles_per_cell;

  cout << "Parameters :" << endl;
  cout << "cell size                  " << cell_size << endl;
  cout << "molecular mass m           " << molecular_mass << endl;
  cout << "number density n           " << n << endl;
  cout << "total # molecules          " << Np << endl;
  cout << "specific weight Fnum       " << fnum << endl;
  cout << "time step from cell size   " << dt0 << endl;
  cout << "time step from d           " << dt1 << endl;
  cout << "time step                  " << dt << endl;
  cout << "mean thermal speed Vn      " << vn << endl;
  cout << "Nstep                      " << nstep << endl;
  cout << "total number of particles  " << num_particles << endl;

#ifdef USE_OCTREE

  //===== OCTREE =====

  // check input parameters!
  OOctree otree(omin,omax,IMAX,JMAX,KMAX,MAX_LEVEL); 

  OREAL cellSize = !otree.realCellSizes(REFINE_MAX);

  cout << "Original octree : num cells " << otree.numCells() << 
    ", num nodes " << otree.numNodes() << ", min cell size " << cellSize << endl;

  // save background cells of octree
  #ifdef SAVE_DEBUG
  {
    TTriangles<OREAL> tris;
    otree.cellsToTriangles(otree.cells(), tris);
    tris.saveSTL(DEBUG_DIR + "originalcells.stl","Octree",true);
  }
  #endif

  // final list of smallest raster cells; this is TEMPORARY list outside
  // octree cells
  OCELLS rastercells;

  // rasterise geometry : refine cells cut by geometry surface
  for (int level = 0; level < REFINE_MAX; level++)
  {
    OCELLS cells;
    otree.rasterizeAtLevel(stator_tris,level,true,cells);

    #ifdef SAVE_DEBUG
    {
      TTriangles<OREAL> tris;
      otree.cellsToTriangles(cells,tris);
      tris.saveSTL(DEBUG_DIR + "rasterised to level " + to_string(level) + ".stl","Octree",true);
    }
    #endif

    // extend rasterised cells by layers of neighbours of same size
    // two calls are necessary to make continuous faces on inner-raster boundary
    otree.addSameSizeNeighbours(cells);
    otree.addSameSizeNeighbours(cells);

    #ifdef SAVE_DEBUG
    {
      TTriangles<OREAL> tris;
      otree.cellsToTriangles(cells,tris);
      tris.saveSTL(DEBUG_DIR + "rasterised extended to level " + to_string(level) + ".stl","Octree",true);
    }
    #endif

    #ifdef _DEBUG
    for (auto c : cells)
    {
      assert(otree.level(c.first) == level);
    }
    #endif

    // refine extended list of raster cells
    otree.refineCells(cells);

    cout << "Refinement at level " << level << ", num cells " << otree.numCells() << 
      ", num nodes " << otree.numNodes() << endl;

    #ifdef SAVE_DEBUG
    {
      TTriangles<OREAL> tris;
      otree.cellsToTriangles(otree.cells(),tris);
      tris.saveSTL(DEBUG_DIR + "octree refined to level " + to_string(level) + ".stl","Octree",true);
    }
    #endif

    if (level == REFINE_MAX - 1)
    {
      cout << "Final rasterisation ..." << endl;
      otree.rasterizeAtLevel(stator_tris,level + 1,true,rastercells);

      #ifdef SAVE_DEBUG
      {
        TTriangles<OREAL> tris;
        otree.cellsToTriangles(otree.cells(),tris);
        tris.saveSTL(DEBUG_DIR + "octree refined to level " + to_string(level + 1) + ".stl","Octree",true);
      }
      #endif

      //#ifdef _DEBUG
      //for (auto c : rastercells)
      //{
      //  assert(c.second.data.second >= 0);
      //}
      //#endif

      #ifdef SAVE_DEBUG
      {
        TTriangles<OREAL> tris;
        otree.cellsToTriangles(rastercells,tris);
        tris.saveSTL(DEBUG_DIR + "raster final at level " + to_string(REFINE_MAX) + ".stl","Octree",true);
      }
      #endif

      // add a level of same size neighbours
      otree.addSameSizeNeighbours(rastercells); 

      #ifdef SAVE_DEBUG
      {
        TTriangles<OREAL> tris;
        otree.cellsToTriangles(rastercells,tris);
        tris.saveSTL(DEBUG_DIR + "raster final extended at level " + to_string(REFINE_MAX) + ".stl","Octree",true);
      }
      #endif

      //#ifdef _DEBUG
      //for (auto c : rastercells)
      //{
      //  assert(c.second.data.second >= 0);
      //}
      //#endif

    } // last level
  } // loop over resterisation levels


  t2 = high_resolution_clock::now();
  duration<double> time_span1 = duration_cast<duration<double>>(t2 - t1);
  cout << "rasterisation and refinement, time " << time_span1.count() << endl; 
  
  // we now process final list of raster cells, classify them by type in its data such as
  // 3 - this is raster cell
  // 1 - inside
  // 2 - outside

  cout << "Classifying cells ..." << endl;
  otree.classifyCells(rastercells);

  t2 = high_resolution_clock::now();
  time_span1 = duration_cast<duration<double>>(t2 - t1);
  cout << "rasterisation, refinement, in/out time " << time_span1.count() << endl; 

  #ifdef SAVE_DEBUG
  {
    TTriangles<OREAL> tris;
    otree.cellsToTriangles(otree.cells(),3,tris);
    tris.saveSTL(DEBUG_DIR + "raster cells.stl","Octree",true);
  }
  #endif

  // get correct in and out indices (outer cells must have bigger size)
  int inindex = 1;
  int outindex = 2;

  // sizes of inner and outer cells
  std::pair<OVECTOR,OVECTOR> minmax1 =  otree.minmax(otree.cells(),1);
  std::pair<OVECTOR,OVECTOR> minmax2 =  otree.minmax(otree.cells(),2);
  OREAL d1 = !(minmax1.second - minmax1.first);
  OREAL d2 = !(minmax2.second - minmax2.first);

  cout << "d1 " << d1 << " d2 " << d2 << endl;

  if (d1 > d2)
  {
    int temp = inindex; inindex = outindex; outindex = temp;
  }

  cout << "Body fitting ..." << endl;

  #ifdef SAVE_DEBUG
  {
    TTriangles<OREAL> tris;
    otree.cellsToTriangles(otree.cells(),-1,tris);
    tris.saveSTL(DEBUG_DIR + "octree cells undefined.stl","Octree",true);
  }
  #endif

  #ifdef SAVE_DEBUG
  {
    TTriangles<OREAL> tris;
    otree.cellsToTriangles(otree.cells(),1,tris);
    tris.saveSTL(DEBUG_DIR + "octree cells type 1.stl","Octree",true);
  }
  #endif

  #ifdef SAVE_DEBUG
  {
    TTriangles<OREAL> tris;
    otree.cellsToTriangles(otree.cells(),2,tris);
    tris.saveSTL(DEBUG_DIR + "octree cells type 2.stl","Octree",true);
  }
  #endif

  #ifdef SAVE_DEBUG
  {
    TTriangles<OREAL> tris;
    otree.cellsToTriangles(otree.cells(),3,tris);
    tris.saveSTL(DEBUG_DIR + "octree cells type 3.stl","Octree",true);
  }
  #endif

  //!!!!!!! signed distance not working well
#if 0
  cout << "Calculating signed distance" << endl;

  // calculate signed siatnce at nodes
  otree.calcSignedDistance();

  cout << "Polygonising surface with marching cubes" << endl;

  // save generated geometry
  {
    TTriangles<OREAL> tris;
    otree.polygoniseDistanceFields(tris,rastercells);
    tris.saveSTL(DEBUG_DIR + BLADE_OUTFILENAME,"Geometry",true);
  }
#endif

#endif // USE_OCTREE

  //===== PARTICLES =====

  std::vector<Particle<T>> particles(num_particles);

  for (int i = 0; i < num_particles; i++)
  {
    initialize_inlet_particle(particles[i],entry_tris,entry_speed_z,
      gen,temperature,m(molar_mass));
  }

  cout << endl << "Press ENTER to exit..." << endl;
  getchar();
}
