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

//#include "Vector.h"
//#include "TMatrix.h"
#include "Conforming3D.h"
#include "OOctree.h"

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

/* INPUT!!!! ******************************************************************/

// save debug info
#define SAVE_DEBUG

// blade, located in the project folder; Z-axis is along span, X - along chord
//#define BLADE_FILENAME "NACA.stl"
#define BLADE_FILENAME "TM pump stator.001.STL"

// octree parameters
#define IMAX 4
#define JMAX 3
#define KMAX 2
#define MAX_LEVEL 10

// raster cells will refined to REFINE_MAX, very important!
#define REFINE_MAX 6

// blade represented by signed distance field
#define BLADE_OUTFILENAME "NACA_SDF.stl"

/******************************************************************************/

using namespace tcad;
using namespace std;
using namespace std::chrono;

int main()
{
  // all about doubles
  #define T double

  // start hight-res timer
  high_resolution_clock::time_point t1,t2;
  t1 = high_resolution_clock::now();

  TTriangles<OREAL> atris;

#if 1

  TBrep<T> rotor,stator;

  bool ok = makeTMPump(0.12,0.05,true, // if false, non-periodic is slow due to multiple boolean operations when adding blades
    0.0007,
    //0.062, // ~= PI20 * R / rotor_numblades
    //0.013, // ~= PI20 * Rshaft / stator_numblades
    -1.0,-1.0, // let it set chords itself
    rotor,stator,
    12,2,
    TMPAirfoil<T>,
    TMPContour<T>,
    TMPAirfoil<T>,
    TMPContour<T>,
    -30.0,+45.0);

  ASSERT(ok);

#else

  // load geometry STL file (solid) 
  std::string partname;
  bool binary;

  // global tolerance
  OREAL tolerance = 0.0;

  bool ok = atris.loadSTL(BLADE_FILENAME,partname,binary,tolerance,OREAL(0.0));
  assert(ok);

  // it is currently size * 0.000001, make 0.00000001
  tolerance *= 0.01;

  #ifdef SAVE_DEBUG
    atris.saveSTL("original.stl","Octree",true);
  #endif

  // exclude repeating nodes
  atris.buildConnectivityArray(tolerance);
#endif

  // calculate cuboid and create octree
  std::pair<OVECTOR,OVECTOR> minmax = atris.minmax();

  OVECTOR min = minmax.first;
  OVECTOR max = minmax.second;
  OVECTOR d = (max - min) * OREAL(0.5);
  OVECTOR c = (max + min) * OREAL(0.5);
  OVECTOR omin = c - d * OREAL(1.2);
  OVECTOR omax = c + d * OREAL(1.2);

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
    tris.saveSTL("originalcells.stl","Octree",true);
  }
  #endif

  // final list of smallest raster cells; this is TEMPORARY list outside
  // octree cells
  OCELLS rastercells;

  // rasterise geometry : refine cells cut by geometry surface
  for (int level = 0; level < REFINE_MAX; level++)
  {
    OCELLS cells;
    otree.rasterizeAtLevel(atris,level,true,cells);

    #ifdef SAVE_DEBUG
    {
      TTriangles<OREAL> tris;
      otree.cellsToTriangles(cells,tris);
      tris.saveSTL("rasterised to level " + to_string(level) + ".stl","Octree",true);
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
      tris.saveSTL("rasterised extended to level " + to_string(level) + ".stl","Octree",true);
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
      tris.saveSTL("octree refined to level " + to_string(level) + ".stl","Octree",true);
    }
    #endif

    if (level == REFINE_MAX - 1)
    {
      cout << "Final rasterisation ..." << endl;
      otree.rasterizeAtLevel(atris,level + 1,true,rastercells);

      #ifdef SAVE_DEBUG
      {
        TTriangles<OREAL> tris;
        otree.cellsToTriangles(otree.cells(),tris);
        tris.saveSTL("octree refined to level " + to_string(level + 1) + ".stl","Octree",true);
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
          tris.saveSTL("raster final at level " + to_string(REFINE_MAX) + ".stl","Octree",true);
        }
      #endif

      // add a level of same size neighbours
      otree.addSameSizeNeighbours(rastercells); 

      #ifdef SAVE_DEBUG
      {
        TTriangles<OREAL> tris;
        otree.cellsToTriangles(rastercells,tris);
        tris.saveSTL("raster final extended at level " + to_string(REFINE_MAX) + ".stl","Octree",true);
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
    tris.saveSTL("raster cells.stl","Octree",true);
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
    tris.saveSTL("octree cells undefined.stl","Octree",true);
  }
  #endif

  #ifdef SAVE_DEBUG
  {
    TTriangles<OREAL> tris;
    otree.cellsToTriangles(otree.cells(),1,tris);
    tris.saveSTL("octree cells type 1.stl","Octree",true);
  }
  #endif

  #ifdef SAVE_DEBUG
  {
    TTriangles<OREAL> tris;
    otree.cellsToTriangles(otree.cells(),2,tris);
    tris.saveSTL("octree cells type 2.stl","Octree",true);
  }
  #endif

  #ifdef SAVE_DEBUG
  {
    TTriangles<OREAL> tris;
    otree.cellsToTriangles(otree.cells(),3,tris);
    tris.saveSTL("octree cells type 3.stl","Octree",true);
  }
  #endif

  cout << "Calculating signed distance" << endl;

  // calculate signed siatnce at nodes
  otree.calcSignedDistance();

  cout << "Polygonising surface with marching cubes" << endl;

  // save generated geometry
  {
    TTriangles<OREAL> tris;
    otree.polygoniseDistanceFields(tris,rastercells);
    tris.saveSTL(BLADE_OUTFILENAME,"Geometry",true);
  }

  //cout << endl << "Press ENTER to exit..." << endl;
  //getchar();
}
