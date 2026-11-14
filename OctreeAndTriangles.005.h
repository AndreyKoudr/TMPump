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

/** Group all triangles for all surface types. object is rotor or stator. */
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
  std::map<IPosition,std::vector<std::pair<int,std::array<TPoint<T>,3>>>,IPosCompare> &cell_tris)
{
  // loop over triangle group types
  for (auto t : object_tris)
  {
    // avoid repetitions, avoid coincidence of
    // cell + tri group type + tri number
    std::set<std::array<LINT,2>> busy;

    // keep tris
    TTriangles<T> &tris = t.second;

    for (LINT j = 0; j < tris.numFaces(); j++)
    {
      std::array<TPoint<T>,3> corners = tris.threeCorners(j);

      std::vector<TPoint<T>> points;

#if 1
      makeTrianglePoints<T>(corners,cell_size,points);
#else
      //TTriangles<T> ttris,tsubtris;
      //makeTrianglePoints<T>(corners,cellSize,points,&ttris,&tsubtris);
      //ttris.saveSTL(DEBUG_DIR + DEBUG_DIR + "tris.stl","Octree",true);
      //tsubtris.saveSTL(DEBUG_DIR + DEBUG_DIR + "subtris.stl","Octree",true);
#endif

      for (auto &p : points)
      {
        IPosition centre;
        if (otree.findCellCentre(p,refine_level,centre)) 
        {
          if (otree.cells()[centre].data.first == cell_type)
          {
            LINT cell_index = otree.positionToBackCellIndex(centre); //!!!!!!!! incorrect
            if (busy.find({cell_index,(LINT) j}) == busy.end())
            {
              cell_tris[centre].push_back({t.first,corners});
              busy.insert({cell_index,(LINT) j});
            }

            otree.cells()[centre].data.first = 4; //!!!!!!!!
          }
        }
      }
    }
  }
}

/** Randomly fill octree cells. */
template <typename T> void fillInnerCells(OOctree &otree, int inindex,
  OOctree &rotree, int inrotorindex, int refine_level,
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
        LIMIT_MIN(n,2);

        // make it random, actual number per cell
        int na = int(T(gen()) / T(gen.max()) * T(n * 2)); // * 2 to get n mean

        // box corners
        std::array<TPoint<T>,8> corners = otree.cell8RealCoordinates(c.first);

        for (int i = 0; i < na; i++)
        {
          // random coordinates
          T E1 = T(gen()) / T(gen.max());
          T E2 = T(gen()) / T(gen.max());
          T E3 = T(gen()) / T(gen.max());
        
          std::array<T,8> func;
          HexLinear3D<T,TPoint<T>>::shapeFunctions(TPoint<T>(E1,E2,E3),func);

          TPoint<T> p;
          for (int i = 0; i < 8; i++)
            p += corners[i] * func[i];

          // check if the particle is inside or in raster of rotor blade
          IPosition rcentre;
          if (rotree.findCellCentre(p,refine_level,rcentre)) 
          {
            if (rotree.cells()[rcentre].data.first != inrotorindex &&
              rotree.cells()[rcentre].data.first != 3)
            {
              Particle<T> particle;
              particle.position = p;
              particle.active = true;
              particles.push_back(particle);
            }
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

  for (auto i = 0; i < particles.size(); i++)
  {
    if (i % decimate == 0)
    {
      TTransform<T> t;
      t.Translate(particles[i].position);

      TTriangles<T> st = stris;
      st.makeTransform(t);

      tris.addTriangles(st);
    }
  }
}
