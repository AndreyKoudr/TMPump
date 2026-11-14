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

#pragma once
     
#include <iostream>

#include "tbasics.h"
#include "ttriangles.h"
#include "tbrep.h"
#include "OOctree.h"
#include "Particle.h"
#include "HexLinear3D.h"

using namespace std;
using namespace tcad;

// save debug info in files
//#define DEBUG_OTRIS

#ifdef NDEBUG
  //!!!!!!! #undef DEBUG_OTRIS
#endif

/** Make octree from triangles. */
template <typename T> void makeOctreeFromTriangles(TTriangles<T> &atris, OOctree &otree, 
  const T cell_size, const int refineLevel, int &inindex, int &outindex,
  const int maxOctreeLevel = 10, T expandCoef = 1.1)
{
  // calculate cuboid and create octree
  std::pair<TPoint<T>,TPoint<T>> minmax = atris.minmax();

  TPoint<T> min = minmax.first;
  TPoint<T> max = minmax.second;
  TPoint<T> dim = (max - min) * 0.5;
  TPoint<T> c = (max + min) * 0.5;
  TPoint<T> omin = c - dim * expandCoef;
  TPoint<T> omax = c + dim * expandCoef;

  // calculate domain size in cells
  TPoint<T> dom = omax - omin;
  int Isize = ROUND(dom.X / cell_size); 
  LIMIT_MIN(Isize,1);
  int Jsize = ROUND(dom.Y / cell_size); 
  LIMIT_MIN(Jsize,1);
  int Ksize = ROUND(dom.Z / cell_size); 
  LIMIT_MIN(Ksize,1);
  int num_cells = Isize * Jsize * Ksize;

  // create background cells (level 0)
  otree.init(omin,omax,Isize,Jsize,Ksize,maxOctreeLevel); 
  otree.createBackground();

  OREAL cellSize = otree.meanRealCellSize(0);

  // save background cells of octree
#ifdef DEBUG_OTRIS
  {
    cout << "Original octree : num cells " << otree.numCells() << 
      ", num nodes " << otree.numNodes() << ", background cell size " << cellSize << endl;

    TTriangles<OREAL> tris;
    otree.cellsToTriangles(otree.cells(), tris);
    tris.saveSTL(DEBUG_DIR + "originalcells.stl","Octree",true);
  }
#endif

  // final list of smallest raster cells; this is TEMPORARY list outside
  // octree cells
  OCELLS rastercells;

  // rasterise geometry : refine cells cut by geometry surface
  for (int level = 0; level < refineLevel; level++)
  {
    OCELLS cells;
    otree.rasterizeAtLevel(atris,level,true,cells);

  #ifdef DEBUG_OTRIS
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

  #ifdef DEBUG_OTRIS
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

  #ifdef DEBUG_OTRIS
    {
      cout << "Refinement at level " << level << ", num cells " << otree.numCells() << 
        ", num nodes " << otree.numNodes() << endl;

      TTriangles<OREAL> tris;
      otree.cellsToTriangles(otree.cells(),tris);
      tris.saveSTL(DEBUG_DIR + "octree refined to level " + to_string(level) + ".stl","Octree",true);
    }
  #endif

    if (level == refineLevel - 1)
    {
  #ifdef DEBUG_OTRIS
      cout << "Final rasterisation ..." << endl;
  #endif
      otree.rasterizeAtLevel(atris,level + 1,true,rastercells);

  #ifdef DEBUG_OTRIS
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

  #ifdef DEBUG_OTRIS
        {
          TTriangles<OREAL> tris;
          otree.cellsToTriangles(rastercells,tris);
          tris.saveSTL(DEBUG_DIR + "raster final at level " + to_string(refineLevel) + ".stl","Octree",true);
        }
  #endif

      // add a level of same size neighbours
      otree.addSameSizeNeighbours(rastercells); 

  #ifdef DEBUG_OTRIS
      {
        TTriangles<OREAL> tris;
        otree.cellsToTriangles(rastercells,tris);
        tris.saveSTL(DEBUG_DIR + "raster final extended at level " + to_string(refineLevel) + ".stl","Octree",true);
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

  // we now process final list of raster cells, classify them by type in its data such as
  // 3 - this is raster cell
  // 1 - inside
  // 2 - outside

#ifdef SAVE_DEBUG
  cout << "Classifying cells ..." << endl;
#endif

  otree.classifyCells(rastercells);

#ifdef DEBUG_OTRIS
  {
    TTriangles<OREAL> tris;
    otree.cellsToTriangles(otree.cells(),3,tris);
    tris.saveSTL(DEBUG_DIR + "raster cells.stl","Octree",true);
  }
#endif

  // get correct in and out indices (outer cells must have bigger size)
  inindex = 1;
  outindex = 2;

  // sizes of inner and outer cells
  std::pair<TPoint<T>,TPoint<T>> minmax1 = otree.minmax(otree.cells(),1);
  std::pair<TPoint<T>,TPoint<T>> minmax2 = otree.minmax(otree.cells(),2);
  OREAL d1 = !(minmax1.second - minmax1.first);
  OREAL d2 = !(minmax2.second - minmax2.first);

  if (d1 > d2)
  {
    int temp = inindex; inindex = outindex; outindex = temp;
  }

#ifdef DEBUG_OTRIS
  {
    TTriangles<OREAL> tris;
    otree.cellsToTrianglesByLevel(otree.cells(),0,tris);
    tris.saveSTL(DEBUG_DIR + "octree background cells.stl","Octree",true);
  }
#endif

#ifdef DEBUG_OTRIS
  {
    TTriangles<OREAL> tris;
    otree.cellsToTriangles(otree.cells(),-1,tris);
    tris.saveSTL(DEBUG_DIR + "octree cells undefined.stl","Octree",true);
  }
#endif

#ifdef DEBUG_OTRIS
  {
    TTriangles<OREAL> tris;
    otree.cellsToTriangles(otree.cells(),inindex,tris);
    tris.saveSTL(DEBUG_DIR + "octree cells inner.stl","Octree",true);
  }
#endif

#ifdef DEBUG_OTRIS
  {
    TTriangles<OREAL> tris;
    otree.cellsToTriangles(otree.cells(),outindex,tris);
    tris.saveSTL(DEBUG_DIR + "octree cells outer.stl","Octree",true);
  }
#endif

#ifdef DEBUG_OTRIS
  {
    TTriangles<OREAL> tris;
    otree.cellsToTriangles(otree.cells(),3,tris);
    tris.saveSTL(DEBUG_DIR + "octree cells raster.stl","Octree",true);
  }
#endif
}

/** Is point inside or outside? Cells of all levels are checked. */
template <typename T> bool pointInsideOutside(TPoint<T> point, OOctree &otree, int inoutindex, 
  int max_refine_level)
{
  for (int i = 0; i <= max_refine_level; i++)
  {
    IPosition centre;
    if (otree.findCellCentre(point,i,centre) && otree.isLeaf(centre) &&
      otree.cells()[centre].data.first == inoutindex)
      return true;
  }
  return false;
}

/** Group all triangles for all surface types. object is stator or stator. */
template <typename T> void groupTriangles(TBrep<T> &object,
  const std::map<int,std::vector<int>> &object_faces, 
  const std::string name,
  const int decimate,
  std::map<int,TTriangles<T>> &object_tris)
{
  // group triangles for every face type
  for (auto t : object_faces)
  {
    TTriangles<T> typetris;
    for (int i = 0; i < int(t.second.size()); i++)
    {
      TTriangles<T> tris;
      object.makeFaceTriangles(t.second[i],tris,decimate,-1.0,0.1,true);
      typetris.addTriangles(tris);
    }

    object_tris[t.first] = typetris;

#ifdef SAVE_DEBUG
    typetris.saveSTL(DEBUG_DIR + DEBUG_DIR + name + "_faces_" + facetypeNames[t.first] + ".stl",name,true);
#endif
  }
}

/** Distribute triangles over octree cells. */
template <typename T> void distributeTrisOverCells(
  const std::map<int,TTriangles<T>> &object_tris,
  OOctree &otree, int refine_level, int cell_type,
  const T cell_size,
                                     // type     face
  std::map<IPosition,std::set<std::pair<int,     int>>,IPosCompare> &cell_tris)
{
  // loop over triangle group types
  for (auto t : object_tris)
  {
    // keep tris
    TTriangles<T> &tris = t.second;

    for (LINT j = 0; j < tris.numFaces(); j++)
    {
      std::array<TPoint<T>,3> corners = tris.threeCorners(j);

      std::vector<TPoint<T>> points;

      // make lots of triangle points
      makeTrianglePoints<T>(corners,cell_size,points);

      // test if any of points is in the cell
      for (auto &p : points)
      {
        IPosition centre;
        if (otree.findCellCentre(p,refine_level,centre)) 
        {
          if (otree.cells()[centre].data.first == cell_type)
          {
            cell_tris[centre].insert({t.first,int(j)});

            // mark this cell as cell number 4 type!
            otree.cells()[centre].data.first = 4; 
          }
        }
      }
    }
  }
}

/** Randomly fill octree cells. */
template <typename T> void fillInnerCells(OOctree &otree, int inindex,
  OOctree &rotree, int outrotorindex, int refine_level,
  std::vector<Particle<T>> &particles, int particles_per_cell, std::mt19937& gen)
{
  for (auto c : otree.cells())
  {
    if (otree.isLeaf(c.first))
    {
      if (c.second.data.first == inindex)
      {
        // 0 - base level, particles_per_cell is specified for it
        auto l = otree.level(c.first);
        // 8 times less for level 1
        int n = particles_per_cell / ((l + 1) * (l + 1) * (l + 1));
        LIMIT_MIN(n,1);
      //!!!!!!!!  LIMIT_MIN(n,2);

        // make it random, actual number per cell
        int na = int(T(gen()) / T(gen.max()) * T(n)); 
   //!!!!!!!!     int na = int(T(gen()) / T(gen.max()) * T(n * 2)); // * 2 to get n mean

        // box corners
        std::array<TPoint<T>,8> corners = otree.cell8RealCoordinates(c.first);

        for (int i = 0; i < na; i++)
        {
          // random coordinates
          T E1 = T(gen()) / T(gen.max()) * 2.0 - 1.0;
          T E2 = T(gen()) / T(gen.max()) * 2.0 - 1.0;
          T E3 = T(gen()) / T(gen.max()) * 2.0 - 1.0;
        
          std::array<T,8> func;
          HexLinear3D<T,TPoint<T>>::shapeFunctions(TPoint<T>(E1,E2,E3),func);

          TPoint<T> p;
          for (int i = 0; i < 8; i++)
            p += corners[i] * func[i];

          // check if the particle is outside rotor blade
          IPosition rcentre;
          if (rotree.insideBox(p)) 
          {
            if (pointInsideOutside(p,rotree,outrotorindex,refine_level))
            {
              Particle<T> particle;
              particle.position = particle.oldposition = p;
              particle.active = true;
              particles.push_back(particle);
            }
          } else
          {
            // outside rotor octree
            Particle<T> particle;
            particle.position = particle.oldposition = p;
            particle.active = true;
            particles.push_back(particle);
          }
        }
      }
    }
  }
}

/** Make triangles as spheres particles for display. */
template <typename T> void makeParticleTriangles(std::vector<Particle<T>> &particles, 
  TTriangles<T> &tris, T R, int decimate = 1)
{
  TBrep<T> sphere;
  sphere.makeSphere(R);
  sphere.closeOuterBoundary();
  TTriangles<T> stris;
  sphere.makeTriangles(stris,100);

  int count = 0;
  for (auto i = 0; i < particles.size(); i++)
  {
    if (!particles[i].active)
      continue;

    if (count % decimate == 0)
    {
      TTransform<T> t;
      t.Translate(particles[i].position);

      TTriangles<T> st = stris;
      st.makeTransform(t);

      tris.addTriangles(st);
    }

    count++;
  }
}

/** Process collision of particle with stator blade in a cell; returns true if processed. */
template <typename T> bool processRotorCell(Particle<T> &particle, IPosition centre,
  std::map<IPosition,std::set<std::pair<int,int>>,IPosCompare> &rotor_cell_tris,
  std::map<int,TTriangles<T>> &rotor_tris,
  T RPM, T blade_temperature, T molecular_mass, std::mt19937& gen)
{
  if (rotor_cell_tris.count(centre))
  {
    for (auto &list : rotor_cell_tris[centre])
    {
      // type
      auto type = list.first;
      // 3 corners
      auto c = rotor_tris[type].threeCorners(list.second);
      T U = 0.0;
      TPoint<T> intr;
      bool backfaced = false;
      if (segTriIntersectMT<T>(particle.oldposition,particle.position,c,U,intr,PARM_TOLERANCE,&backfaced))
      {
        if (backfaced)
        {
          // this must not be, exclude this particle, the particle is inside wing
          // what must not happen
          particle.active = false;
        } else
        {
          // out of chamber through exit or entry
          if (type == WING)
          {
            // reflect moving blade
            T r = sqrt(intr.X * intr.X + intr.Y * intr.Y);
            T speed = bladeSpeed<T>(r,RPM);
            T angle = atan2(intr.Y,intr.X);
            TPoint<T> wing_velocity = +(TPoint<T>(0.0,0.0,1.0) ^ TPoint<T>(intr.X,intr.Y)) * speed;

            TPoint<T> normal = +((c[1] - c[0]) ^ (c[2] - c[0]));
            applyParticleDiffuseReflection<T>(particle, 
              wing_velocity,normal,blade_temperature,molecular_mass,gen);
          }
        }
        // processed
        return true;
      }
    }
    return false;
  } else
  {
    return false;
  }
}

/** Process collision of particle with stator blade; returns true if processed. */
template <typename T> bool processRotorBlade(Particle<T> &particle, 
  OOctree &rotree, int refine_level, 
  std::map<IPosition,std::set<std::pair<int,int>>,IPosCompare> &rotor_cell_tris,
  std::map<int,TTriangles<T>> &rotor_tris,
  T RPM, T blade_temperature, T molecular_mass, std::mt19937& gen)
{
  IPosition centre,oldcentre;
  if (rotree.findCellCentre(particle.position,refine_level,centre) &&
    rotree.findCellCentre(particle.oldposition,refine_level,oldcentre))
  {
    if (rotree.cells()[centre].data.first == 4)
    {
      if (processRotorCell<T>(particle,centre,rotor_cell_tris,rotor_tris,
        RPM,blade_temperature,molecular_mass,gen))
        return true;
    } 

    if (rotree.cells()[oldcentre].data.first == 4)
    {
      if (processRotorCell<T>(particle,oldcentre,rotor_cell_tris,rotor_tris,
        RPM,blade_temperature,molecular_mass,gen))
        return true;
    }

    return false;
  } else
  {
    // not processed
    return false;
  }
}

/** Process collision of particle with stator blade in a cell; returns true if processed. */
template <typename T> bool processStatorCell(Particle<T> &particle, IPosition centre,
  std::map<IPosition,std::set<std::pair<int,int>>,IPosCompare> &stator_cell_tris,
  std::map<int,TTriangles<T>> &stator_tris,
  T RPM, T wall_temperature, T molecular_mass, std::mt19937& gen, T period)
{
  if (stator_cell_tris.count(centre))
  {
    for (auto &list : stator_cell_tris[centre])
    {
      // type
      auto type = list.first;
      // corners
      auto c = stator_tris[type].threeCorners(list.second);
      T U = 0.0;
      TPoint<T> intr;
      bool backfaced = false;
      if (segTriIntersectMT<T>(particle.oldposition,particle.position,c,U,intr,PARM_TOLERANCE,&backfaced))
      {
        if (backfaced)
        {
          // this must not be, exclude this particle
          particle.active = false;
        } else
        {
          // out of chamber through exit or entry
          if (type == CYLINDER_TOP || type == CYLINDER_BOTTOM)
          {
            // escape from entry/exit
            particle.active = false;
          } else if (type == WING || type == CYLINDER || type == SHAFT)
          {
            // reflect from stationary wall (all stator surfaces)
            TPoint<T> wall_velocity;
            TPoint<T> normal = +((c[1] - c[0]) ^ (c[2] - c[0]));
            applyParticleDiffuseReflection<T>(particle, 
              wall_velocity,normal,wall_temperature,molecular_mass,gen);
          } else if (type == CYLINDER_FLANK0)
          {
            // jump through period to another periodic wall, rotate CCW around Z
            TTransform<T> t;
            t.Rotate(TPoint<T>(0.0,0.0,1.0),period * CPI);
            particle.position = particle.oldposition = t.applyTransform(particle.position);
            particle.velocity = t.applyTransform(particle.velocity);
          } else if (type == CYLINDER_FLANK1)
          {
            // jump through period to another periodic wall, rotate CW around Z
            TTransform<T> t;
            t.Rotate(TPoint<T>(0.0,0.0,1.0),-period * CPI);
            particle.position = particle.oldposition = t.applyTransform(particle.position);
            particle.velocity = t.applyTransform(particle.velocity);
          }
        }
        return true;
      }
    }
    return false;
  } else
  {
    return false;
  }
}

/** Process collision of particle with stator blade and walls; returns true if processed. */
template <typename T> bool processStator(Particle<T> &particle, 
  OOctree &otree, int refine_level, 
  std::map<IPosition,std::set<std::pair<int,int>>,IPosCompare> &stator_cell_tris,
  std::map<int,TTriangles<T>> &stator_tris,
  T RPM, T wall_temperature, T molecular_mass, std::mt19937& gen, T period)
{
  IPosition centre,oldcentre;
  if (otree.findCellCentre(particle.position,refine_level,centre) &&
    otree.findCellCentre(particle.oldposition,refine_level,oldcentre))
  {
    if (otree.cells()[centre].data.first == 4)
    {
      if (processStatorCell<T>(particle,centre,stator_cell_tris,stator_tris,
        RPM,wall_temperature,molecular_mass,gen,period))
        return true;
    } 

    if (otree.cells()[oldcentre].data.first == 4)
    {
      if (processStatorCell<T>(particle,oldcentre,stator_cell_tris,stator_tris,
        RPM,wall_temperature,molecular_mass,gen,period))
        return true;
    }

    return false;
  } else
  {
    // not processed
    return false;
  }
}

/** Calculate number of particles in cells. */
template <typename T> void countParticlesInCells(std::vector<Particle<T>> &particles, OOctree &otree, 
  std::map<IPosition,int,IPosCompare> &counts, int refine_level = 0)
{
  counts.clear();

  // fill all cells with zero values
  for (auto c : otree.cells())
  {
    if (otree.level(c.first) == refine_level)
    {
      counts[c.first] = 0;
    }
  }

  // now take counts from particles
  for (auto &p : particles)
  {
    if (p.active)
    {
      IPosition centre;
      if (otree.findCellCentre(p.position,refine_level,centre)) 
      {
        counts[centre]++;
      }
    }
  }
}
  
#if 0

/** Export density in cells in VTK. */
template <typename T> bool exportDensityVTK(std::vector<Particle<T>> &particles, OOctree &otree,
  const T Fnum, const T temperature, const std::string &filename)
{
  std::map<IPosition,int,IPosCompare> counts;
  countParticlesInCells<T>(particles,otree,counts,0);

  std::vector<std::string> lines;

  // header
  lines.push_back("# vtk DataFile Version 3.0\n");
  lines.push_back("Point cloud with pressure data\n");
  lines.push_back("ASCII\n");

  lines.push_back("DATASET POLYDATA\n");

  // # points
  auto n = otree.IJKnumCells.X * otree.IJKnumCells.Y * otree.IJKnumCells.Z;

  lines.push_back("POINTS " + to_string(n) + " float\n");

  for (auto c : counts)
  {
    // point coordinate
    TPoint<T> p = otree.intToRealCoord(c.first);
    lines.push_back(to_string(p.X,12) + " " + to_string(p.Y,12) + " " + to_string(p.Z,12) + "\n");
  }

  lines.push_back("\n");

  lines.push_back("VERTICES " + to_string(n) + " " + to_string(n * 2) + "\n");
  for (auto i = 0; i < n; i++)
  {
    lines.push_back("1 "+ to_string(i) + "\n");
  }

  lines.push_back("\n");

  lines.push_back("POINT_DATA " + to_string(n) + "\n");
  lines.push_back("SCALARS pressure float 1\n");
  lines.push_back("LOOKUP_TABLE default\n");

  for (auto c : counts)
  {
    // calculate pressure
    T vol = otree.cellVolume(c.first);
    int n = c.second;
    T p = n * Fnum * KB * temperature / vol; 

    lines.push_back(to_string(p,12) + "\n");
  }

  return writeLines(lines,filename);
}

#else

/** Export density in cells in VTK. */
template <typename T> bool exportDensityVTK(std::vector<Particle<T>> &particles, OOctree &otree,
  const T Fnum, const T temperature, const std::string &filename)
{
  std::map<IPosition,int,IPosCompare> temp;
  countParticlesInCells<T>(particles,otree,temp,0);

  std::map<IPosition,int,IPosCompareZ> counts;
  counts.insert(temp.begin(),temp.end());

  std::vector<std::string> lines;

  // header
  lines.push_back("# vtk DataFile Version 3.0\n");
  lines.push_back("3D Pressure Distribution\n");
  lines.push_back("ASCII\n");
  lines.push_back("DATASET STRUCTURED_POINTS\n");

  // X is the first
  lines.push_back("DIMENSIONS " + to_string(otree.IJKnumCells.X) + " " + 
    to_string(otree.IJKnumCells.Y) + " " + to_string(otree.IJKnumCells.Z) + "\n");

  TPoint<T> cellsize = otree.realCellSizes(0);
  lines.push_back("ORIGIN " + to_string(otree.boxMin.X) + " " + to_string(otree.boxMin.Y) + " " + to_string(otree.boxMin.Z) + "\n");
  lines.push_back("SPACING " + to_string(cellsize.X) + " " + to_string(cellsize.Y) + " " + to_string(cellsize.Z) + "\n");

  lines.push_back("\n");

  lines.push_back("POINT_DATA " + to_string(otree.IJKnumCells.X * 
    otree.IJKnumCells.Y * otree.IJKnumCells.Z) + "\n");

  lines.push_back("SCALARS pressure float\n");
  lines.push_back("LOOKUP_TABLE default\n");

  for (auto c : counts)
  {
    // calculate pressure
    T vol = otree.cellVolume(c.first);
    int n = c.second;
    T p = n * Fnum * KB * temperature / vol; 

    lines.push_back(to_string(p,12) + "\n");
  }

  return writeLines(lines,filename);
}

#endif
  
/** Export density in cells in VTK. */
template <typename T> bool exportDensityCSV(std::vector<Particle<T>> &particles, OOctree &otree,
  const T Fnum, const T temperature, const std::string &filename)
{
  std::map<IPosition,int,IPosCompare> counts;
  countParticlesInCells<T>(particles,otree,counts,0);

  std::vector<std::string> lines;

  // header
  lines.push_back("X,Y,Z,Pressure\n");

  for (auto c : counts)
  {
    TPoint<T> centre = otree.intToRealCoord(c.first);

    // calculate pressure
    T vol = otree.cellVolume(c.first);
    int n = c.second;
    T p = n * Fnum * KB * temperature / vol;

    lines.push_back(to_string(centre.X,12) + "," + 
      to_string(centre.Y,12) + "," + to_string(centre.Z,12) + "," + to_string(p,12) + "\n");
  }

  return writeLines(lines,filename);
}
