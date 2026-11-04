#pragma once

#include <limits>
#include <algorithm>
#include <functional>
#include <numeric> 
#include <array>
#include <string>
#include <vector>
#include <assert.h>
#include "Types.h"
#include "Matrix.h"
#include "GaussInt.h"

/**
  Quadratic finite element
  ------------------------

  J.J.Connor, C.A.Brebbia. Finite Element Techniques for Fluid Flow.

  Local coordinates! : from -1 to +1


         (7)----------18-----------(6)
         /|                        /|    
        / |                       / |    
      19  |                      17 |     
      /   |                     /   |          
     /   15                    /    14         
   (4)-----------16----------(5)    |    
    |     |     W ^ V         |     |     
    |     |       |/          |     |    
    |     |       O-> U       |     |
    |    (3)----------10------|----(2)
   12    /                    13   / 
    |   /                     |   / 
    |  11                     |  9
    | /                       | /
    |/                        |/
   (0)-----------8-----------(1)
       
  T - class for real numbers
  Tvector - class for 3D vector with up to 4 components
*/

template <typename T, typename Tvector> class HexQuadratic3D {
public:
  /** Total number of nodes */
  static const size_t NUM_NODES = 20;

  /** Number of basic nodes which always exist */
  static const size_t NUM_BASICNODES = 8;

  /** Constructor */
  HexQuadratic3D() = delete;

  /** Constructors */
  HexQuadratic3D(const std::array<LINT,NUM_NODES> &nodes);

  /** Copy constructor */
  HexQuadratic3D(const HexQuadratic3D& copy);

  /** Assignment operator */
  HexQuadratic3D& operator=(const HexQuadratic3D& copy);

  /** Nodes */
  std::array<LINT,NUM_NODES>& nodes()
  {
    return nodes_;
  }

  /** Get min and max node numbers */
  void  nodeMinMax(LINT& min, LINT& max);

  /** Compute shape functions, local coordinates (parms) being [-1..+1]! */
  static void shapeFunctions(const Tvector &parms, std::array<T,NUM_NODES> &func);

  /** Compute derivatives of shape functions, local coordinates (parms) [-1..+1]! */
  static void shapeDerivatives(const Tvector &parms, std::array<Tvector,NUM_NODES> &func);

  /** Compute linear shape functions from 8 basic nodes, local coordinates (parms) being [-1..+1]! */
  static void shapeLinearFunctions(const Tvector &parms, std::array<T,NUM_NODES> &func);

  /** Compute derivatives of linear conforming shape functions, local coordinates (parms) [-1..+1],
    calculation is made by calling shapeFunctions() to compute finite differences;
    result is approximate. */
  static void shapeDerivativesSimple(const Tvector &parms, std::array<Tvector,NUM_NODES> &func);

  /** Get parameters of all nodes, -1..+1 */
  static const std::array<Tvector,NUM_NODES>& nodeParms() 
  {
    return nodeParms_;
  }

  /** Element 6 faces, see the picture for elemenFaces_. */
  static const std::array<std::array<int,8>,6>& faces();

  /** Element edges, the third number is middle node [12][3]. */
  static const std::array<std::array<int,3>,12>& edges();

  /** Two faces for every edge [12][2]. */
  static const std::array<std::array<int,2>,12>& edgeFaces();

  /** Return last error string. */
  static std::string errorString()
  {
    return errorString_;
  }

  // Finite element matrices

  /** Get 20 x 3 node coordinate matrix. */
  TMatrix<T,NUM_NODES,3> coordMatrix(const std::vector<Tvector> &totalCoordinates)
  {
    TMatrix<T,NUM_NODES,3> coords;

    for (int i = 0; i < NUM_NODES; i++)
    {
      if (nodes_[i] >= 0)
      {
        coords[i][0] = totalCoordinates[nodes_[i]].XYZ[0];
        coords[i][1] = totalCoordinates[nodes_[i]].XYZ[1];
        coords[i][2] = totalCoordinates[nodes_[i]].XYZ[2];
      } else
      {
        coords[i][0] = coords[i][1] = coords[i][2] = T(0.0);
      }
    }

    return coords;
  }

  /** Get 20 x 1 shape function matrix, zeroes for non-existing nodes. */
  static TMatrix<T,NUM_NODES,1> shapeFunctionMatrix(const std::array<T,NUM_NODES> &func)
  {
    TMatrix<T,NUM_NODES,1> funcmatrix;
    funcmatrix.fill(T0);

    for (int i = 0; i < NUM_NODES; i++)
    {
      funcmatrix[i][0] = func[i];
    }

    return funcmatrix;
  }

  /** Get 20 x 3 shape function derivative matrix. */
  static TMatrix<T,NUM_NODES,3> shapeFunctionDerivativesMatrix(const std::array<Tvector,NUM_NODES> &func)
  {
    TMatrix<T,NUM_NODES,3> dermatrix;
    dermatrix.fill(T0);

    for (int i = 0; i < NUM_NODES; i++)
    {
      dermatrix[i][0] = func[i].XYZ[0];
      dermatrix[i][1] = func[i].XYZ[1];
      dermatrix[i][2] = func[i].XYZ[2];
    }

    return dermatrix;
  }


  /**
    Calculate matrix B to get deriveatives on X,Y,Z (real coordinates) by
    {Ux,Uy,Uz} = B[3,NUM_NODES] * { U nodal values [NUM_NODES] }
    See J.J.Connor, C.A.Brebbia Finite Element Techniques for Fluid Flow, page 125, equation (i).
  */
  bool Bmatrix(const std::vector<Tvector> &totalCoordinates, 
    const std::array<Tvector,NUM_NODES> &derivatives, TMatrix<T,3,NUM_NODES> &B, T &jacobian);

  /** Calculate element stiffness matrix [K]. */
  bool  stiffnessMatrix(const std::vector<Tvector>& totalCoordinates,
    TMatrix<T,HexQuadratic3D<T,Tvector>::NUM_NODES,HexQuadratic3D<T,Tvector>::NUM_NODES>& matrix, const int power);

  /** Calculate element mass matrix [M]. */
  bool  massMatrix(const std::vector<Tvector>& totalCoordinates,
    TMatrix<T,HexQuadratic3D<T,Tvector>::NUM_NODES,HexQuadratic3D<T,Tvector>::NUM_NODES>& matrix, const int power);

  /** Returns volume, -1.0 in case of failure. */
  T volume(const std::vector<Tvector>& totalCoordinates, const int power);

  // ===== Solution vector =====

  /** Get nodal values for this element */
  TMatrix<T,1,HexQuadratic3D<T,Tvector>::NUM_NODES> values(const std::vector<T>& totalValues)
  {
    TMatrix<T,1,HexQuadratic3D<T,Tvector>::NUM_NODES> vs;
    vs.fill(T0);

    for (size_t i = 0; i < HexQuadratic3D<T,Tvector>::NUM_NODES; i++)
    {
      vs[0][i] = totalValues[nodes_[i]];
    }

    return vs;
  }


  /** Calculate coordinate at a point defined by parameter */
  Tvector coordinate(const Tvector& parms, const std::vector<Tvector>& totalCoordinates);

  /** Calculate a value at a point defined by parameter */
  T value(const Tvector& parms, const std::vector<T>& totalValues);

  /** Calculate value gradient at a point defined by parameter */
  Tvector gradient(const Tvector& parms, const std::vector<Tvector> &totalCoordinates,
    const std::vector<T>& totalValues);

private:
  /** Node numbers. */
  std::array<LINT,NUM_NODES> nodes_;

  /** Parametric values for all nodes plus central node */
  static const std::array<Tvector,NUM_NODES> nodeParms_;

  /** Last error string */
  static std::string errorString_;

};


template <typename  T, typename Tvector>
  std::string HexQuadratic3D<T,Tvector>::errorString_;

template <typename  T, typename Tvector> HexQuadratic3D<T,Tvector>::
  HexQuadratic3D(const std::array<LINT,NUM_NODES> &nodes)
{
  nodes_ = nodes;
}

template <typename  T, typename Tvector> HexQuadratic3D<T,Tvector>::
  HexQuadratic3D(const HexQuadratic3D& copy)
{
  nodes_ = copy.nodes_;
}

template <typename  T, typename Tvector> HexQuadratic3D<T,Tvector>& HexQuadratic3D<T,Tvector>::
  operator = (const HexQuadratic3D& copy)
{
  nodes_ = copy.nodes_;

  return *this;
}


template <typename  T, typename Tvector> const std::array<std::array<int,8>,6>& 
  HexQuadratic3D<T,Tvector>::faces()
{

/**
  Element faces :
 
        7-------------------6 
       /|                  /|    
      / |      (5)        / |     
     /  |                /  |          
    4-----------(3)-----5   |    
    |   |               |   |     
    |(0)|               |(1)|    
    |   |               |   |
    |   3----(2)--------|---2 
    |  /                |  / 
    | /        (4)      | /
    |/                  |/
    0-------------------1

  Hex faces nodal numbers, mind the order : 
 -e direction (-x)
 +e direction (+x)
 -n direction (-y)
 +n direction (+y)
 -s direction (-z)
 +s direction (+z)
*/

  static const std::array<std::array<int,8>,6>  elementFaces_ =
  {{
    {3,0,4,7,   11,12,19,15},
    {1,2,6,5,   9,14,17,13},
    {0,1,5,4,   8,13,16,12},
    {2,3,7,6,   10,15,18,14},
    {0,3,2,1,   11,10,9,8},
    {7,4,5,6,   19,16,17,18}
  }};

  return elementFaces_;
}

template <typename  T, typename Tvector> const std::array<std::array<int,3>,12>& 
  HexQuadratic3D<T,Tvector>::edges()
{

/**
  Element edges, the third number is middle edge node [12][3] :

        7--------(6)--------6 
       /|                  /|    
     (7)|                (5)|     
     /  |                /  |          
    4---------(4)-------5  (10)  
    |  (11)             |   |     
    |   |               |   |    
    |   |               |   |
   (8)  3-------(2)----(9)--2 
    |  /                |  / 
    |(3)                |(1) 
    |/                  |/
    0--------(0)--------1
                                          
*/

  static const std::array<std::array<int,3>,12> elementEdges_ =
  {{
    {0,1,8},
    {1,2,9},
    {2,3,10},
    {3,0,11},

    {4,5,16},
    {5,6,17},
    {6,7,18},
    {7,4,19},

    {0,4,12},
    {1,5,13},
    {2,6,14},
    {3,7,15}
  }};

  return elementEdges_;
}

template <typename  T, typename Tvector> const std::array<std::array<int,2>,12>& 
  HexQuadratic3D<T,Tvector>::edgeFaces()
{


/**
  Two faces for every edge.

  Edges :
        7--------(6)--------6 
       /|                  /|    
     (7)|                (5)|     
     /  |                /  |          
    4---------(4)-------5  (10)  
    |  (11)             |   |     
    |   |               |   |    
    |   |               |   |
   (8)  3-------(2)----(9)--2 
    |  /                |  / 
    |(3)                |(1) 
    |/                  |/
    0--------(0)--------1


  Faces :
        7-------------------6 
       /|                  /|    
      / |      (5)        / |     
     /  |                /  |          
    4-----------(3)-----5   |    
    |   |               |   |     
    |(0)|               |(1)|    
    |   |               |   |
    |   3----(2)--------|---2 
    |  /                |  / 
    | /        (4)      | /
    |/                  |/
    0-------------------1

*/

  static const std::array<std::array<int,2>,12> edgeFaces_ =
  {{
    {2,4},
    {1,4},
    {3,4},
    {0,4},

    {2,5},
    {1,5},
    {3,5},
    {0,5},

    {0,2},
    {2,1},
    {1,3},
    {3,0}
  }};

  return edgeFaces_;
}

template <typename  T, typename Tvector> void HexQuadratic3D<T,Tvector>::shapeFunctions(
  const Tvector &parms, std::array<T,NUM_NODES> &func)
{
#ifndef NDEBUG
#if 1
  T tol = TOLERANCE(T);
  assert(parms.XYZ[0] >= -T1 - tol && parms.XYZ[0] <= T1 + tol);
  assert(parms.XYZ[1] >= -T1 - tol && parms.XYZ[1] <= T1 + tol);
  assert(parms.XYZ[2] >= -T1 - tol && parms.XYZ[2] <= T1 + tol);
#else
  T tol = TOLERANCE(T);
  assert(parms.XYZ[0] >= T0 - tol && parms.XYZ[0] <= T1 + tol);
  assert(parms.XYZ[1] >= T0 - tol && parms.XYZ[1] <= T1 + tol);
  assert(parms.XYZ[2] >= T0 - tol && parms.XYZ[2] <= T1 + tol);
#endif
#endif

  std::fill(func.begin(),func.end(),T0);

  T U2 = parms.X * parms.X;
  T V2 = parms.Y * parms.Y;
  T W2 = parms.Z * parms.Z;

  for (int i = 0; i < 8; i++)
  {
    T UUi = parms.X * nodeParms_[i].X; 
    T VVi = parms.Y * nodeParms_[i].Y; 
    T WWi = parms.Z * nodeParms_[i].Z; 

    T part1 = (T1 + UUi) * (T1 + VVi) * (T1 + WWi);
    T part2 = (UUi + VVi + WWi - T2);

    func[i] = part1 * part2 * 0.125;
  }

/*
         (7)----------18-----------(6)
         /|                        /|    
        / |                       / |    
      19  |                      17 |     
      /   |                     /   |          
     /   15                    /    14         
   (4)-----------16----------(5)    |    
    |     |     W ^ V         |     |     
    |     |       |/          |     |    
    |     |       O-> U       |     |
    |    (3)----------10------|----(2)
   12    /                    13   / 
    |   /                     |   / 
    |  11                     |  9
    | /                       | /
    |/                        |/
   (0)-----------8-----------(1)
*/


  static const std::array<int,4> group0 = {8,10,16,18};
  static const std::array<int,4> group1 = {9,11,17,19};
  static const std::array<int,4> group2 = {12,13,14,15};

  for (int i = 0; i < 4; i++)
  {
    func[group0[i]] = 0.25 * (T1 - U2) * 
      (T1 + parms.Y * nodeParms_[group0[i]].Y) * 
      (T1 + parms.Z * nodeParms_[group0[i]].Z);
    func[group1[i]] = 0.25 * (T1 - V2) * 
      (T1 + parms.X * nodeParms_[group1[i]].X) * 
      (T1 + parms.Z * nodeParms_[group1[i]].Z);
    func[group2[i]] = 0.25 * (T1 - W2) * 
      (T1 + parms.X * nodeParms_[group2[i]].X) * 
      (T1 + parms.Y * nodeParms_[group2[i]].Y);
  }
}

template <typename  T, typename Tvector> void HexQuadratic3D<T,Tvector>::shapeLinearFunctions(
  const Tvector &parms, std::array<T,NUM_NODES> &func)
{
#ifndef NDEBUG
#if 1
  T tol = TOLERANCE(T);
  assert(parms.XYZ[0] >= -T1 - tol && parms.XYZ[0] <= T1 + tol);
  assert(parms.XYZ[1] >= -T1 - tol && parms.XYZ[1] <= T1 + tol);
  assert(parms.XYZ[2] >= -T1 - tol && parms.XYZ[2] <= T1 + tol);
#else
  T tol = TOLERANCE(T);
  assert(parms.XYZ[0] >= T0 - tol && parms.XYZ[0] <= T1 + tol);
  assert(parms.XYZ[1] >= T0 - tol && parms.XYZ[1] <= T1 + tol);
  assert(parms.XYZ[2] >= T0 - tol && parms.XYZ[2] <= T1 + tol);
#endif
#endif

  std::fill(func.begin(),func.end(),T0);

  for (int i = 0; i < 8; i++)
  {
    T UUi = parms.X * nodeParms_[i].X; 
    T VVi = parms.Y * nodeParms_[i].Y; 
    T WWi = parms.Z * nodeParms_[i].Z; 

    func[i] = (T1 + UUi) * (T1 + VVi) * (T1 + WWi) * 0.125;
  }
}

template <class T, class Tvector> void HexQuadratic3D<T,Tvector>::shapeDerivatives(
  const Tvector &parms, std::array<Tvector,NUM_NODES> &func)
{
#ifndef NDEBUG
#if 1
  T tol = TOLERANCE(T);
  assert(parms.XYZ[0] >= -T1 - tol && parms.XYZ[0] <= T1 + tol);
  assert(parms.XYZ[1] >= -T1 - tol && parms.XYZ[1] <= T1 + tol);
  assert(parms.XYZ[2] >= -T1 - tol && parms.XYZ[2] <= T1 + tol);
#else
  T tol = TOLERANCE(T);
  assert(parms.XYZ[0] >= T0 - tol && parms.XYZ[0] <= T1 + tol);
  assert(parms.XYZ[1] >= T0 - tol && parms.XYZ[1] <= T1 + tol);
  assert(parms.XYZ[2] >= T0 - tol && parms.XYZ[2] <= T1 + tol);
#endif
#endif

  std::fill(func.begin(),func.end(),T0);

  T U2 = parms.X * parms.X;
  T V2 = parms.Y * parms.Y;
  T W2 = parms.Z * parms.Z;

  for (int i = 0; i < 8; i++)
  {
    T UUi = parms.X * nodeParms_[i].X; 
    T VVi = parms.Y * nodeParms_[i].Y; 
    T WWi = parms.Z * nodeParms_[i].Z; 

    // part1(U,V,W) ((T1 + UUi) * (T1 + VVi) * (T1 + WWi))
    T part1 = (T1 + UUi) * (T1 + VVi) * (T1 + WWi);

    Tvector part1d; 
    part1d.X = nodeParms_[i].X * (T1 + VVi) * (T1 + WWi);
    part1d.Y = (T1 + UUi) * nodeParms_[i].Y * (T1 + WWi);
    part1d.Z = (T1 + UUi) * (T1 + VVi) * nodeParms_[i].Z;

    // part2(U,V,W) (UUI + VVi + WWi - T2)
    T part2 = (UUi + VVi + WWi - T2);

    Tvector part2d; 
    part2d.X = nodeParms_[i].X;
    part2d.Y = nodeParms_[i].Y;
    part2d.Z = nodeParms_[i].Z;

    func[i] = (part1d * part2 + part2d * part1) * 0.125;
  }

/*
         (7)----------18-----------(6)
         /|                        /|    
        / |                       / |    
      19  |                      17 |     
      /   |                     /   |          
     /   15                    /    14         
   (4)-----------16----------(5)    |    
    |     |     W ^ V         |     |     
    |     |       |/          |     |    
    |     |       O-> U       |     |
    |    (3)----------10------|----(2)
   12    /                    13   / 
    |   /                     |   / 
    |  11                     |  9
    | /                       | /
    |/                        |/
   (0)-----------8-----------(1)
*/


  static const std::array<int,4> group0 = {8,10,16,18};
  static const std::array<int,4> group1 = {9,11,17,19};
  static const std::array<int,4> group2 = {12,13,14,15};

  for (int i = 0; i < 4; i++)
  {
    //func[group0[i]] = 0.25 * (T1 - parms.X * parms.X) * 
    //  (T1 + parms.Y * nodeParms_[group0[i]].Y) * 
    //  (T1 + parms.Z * nodeParms_[group0[i]].Z);
    int index0 = group0[i];
    func[index0].X = 0.25 * (-2.0 * parms.X) * 
      (T1 + parms.Y * nodeParms_[index0].Y) * 
      (T1 + parms.Z * nodeParms_[index0].Z);
    func[index0].Y = 0.25 * (T1 - parms.X * parms.X) * 
      (nodeParms_[index0].Y) * 
      (T1 + parms.Z * nodeParms_[index0].Z);
    func[index0].Z = 0.25 * (T1 - parms.X * parms.X) * 
      (T1 + parms.Y * nodeParms_[index0].Y) * 
      (nodeParms_[index0].Z);

    //func[group1[i]] = 0.25 * (T1 - parms.Y * parms.Y) * 
    //  (T1 + parms.X * nodeParms_[group1[i]].X) * 
    //  (T1 + parms.Z * nodeParms_[group1[i]].Z);
    int index1 = group1[i];
    func[index1].X = 0.25 * (T1 - parms.Y * parms.Y) * 
      (nodeParms_[index1].X) * 
      (T1 + parms.Z * nodeParms_[index1].Z);
    func[index1].Y = 0.25 * (-2.0 * parms.Y) * 
      (T1 + parms.X * nodeParms_[index1].X) * 
      (T1 + parms.Z * nodeParms_[index1].Z);
    func[index1].Z = 0.25 * (T1 - parms.Y * parms.Y) * 
      (T1 + parms.X * nodeParms_[index1].X) * 
      (nodeParms_[index1].Z);

    //func[group2[i]] = 0.25 * (T1 - parms.Z * parms.Z) * 
    //  (T1 + parms.X * nodeParms_[group2[i]].X) * 
    //  (T1 + parms.Y * nodeParms_[group2[i]].Y);
    int index2 = group2[i];
    func[index2].X = 0.25 * (T1 - parms.Z * parms.Z) * 
      (nodeParms_[index2].X) * 
      (T1 + parms.Y * nodeParms_[index2].Y);
    func[index2].Y = 0.25 * (T1 - parms.Z * parms.Z) * 
      (T1 + parms.X * nodeParms_[index2].X) * 
      (nodeParms_[index2].Y);
    func[index2].Z = 0.25 * (-2.0 * parms.Z) * 
      (T1 + parms.X * nodeParms_[index2].X) * 
      (T1 + parms.Y * nodeParms_[index2].Y);
  }
}

template <class T, class Tvector> void HexQuadratic3D<T,Tvector>::shapeDerivativesSimple(
  const Tvector &parms, std::array<Tvector,NUM_NODES> &func)
{
  T step = T(0.001);
  Tvector coords0 = parms - Tvector(step,step,step);
  Tvector coords1 = parms + Tvector(step,step,step);

  Tvector dcoords;
  for (int j = 0; j < 3; j++)
  {
#if 1
    if (coords0.XYZ[j] < -1.0)
      coords0.XYZ[j] = -1.0;
    if (coords1.XYZ[j] > 1.0)
      coords1.XYZ[j] = 1.0;
#else
    if (coords0.XYZ[j] < 0)
      coords0.XYZ[j] = 0;
    if (coords1.XYZ[j] > 1)
      coords1.XYZ[j] = 1;
#endif
    dcoords.XYZ[j] = coords1.XYZ[j] - coords0.XYZ[j];
    assert(dcoords.XYZ[j] > 0);
  }

  for (int j = 0; j < 3; j++)
  {
    std::array<T,NUM_NODES> func0;
    Tvector c0 = parms;
    c0.XYZ[j] = coords0[j];

    shapeFunctions(c0,func0);

    std::array<T,NUM_NODES> func1;
    Tvector c1 = parms;
    c1.XYZ[j] = coords1[j];
        
    shapeFunctions(c1,func1);

    for (int k = 0; k < NUM_NODES; k++)
    {
      func[k].XYZ[j] = (func1[k] - func0[k]) / dcoords[j];
    }
  }
}

template <typename T, typename Tvector> const std::array<Tvector,20> 
  HexQuadratic3D<T,Tvector>::nodeParms_ =
{
  Tvector(-1,-1,-1),
  Tvector(+1,-1,-1),
  Tvector(+1,+1,-1),
  Tvector(-1,+1,-1),

  Tvector(-1,-1,+1),
  Tvector(+1,-1,+1),
  Tvector(+1,+1,+1),
  Tvector(-1,+1,+1),

  Tvector( 0,-1,-1),
  Tvector(+1, 0,-1),
  Tvector( 0,+1,-1),
  Tvector(-1, 0,-1),

  Tvector(-1,-1, 0),
  Tvector(+1,-1, 0),
  Tvector(+1,+1, 0),
  Tvector(-1,+1, 0),

  Tvector( 0,-1,+1),
  Tvector(+1, 0,+1),
  Tvector( 0,+1,+1),
  Tvector(-1, 0,+1)
};

// Multiply shape functions by values at nodes and sum up
template <typename T, typename Tvalue, size_t N> Tvalue shapeFuncValue(const std::array<T,N> &func,
  const std::array<Tvalue,N> &values, Tvalue zero)
{
  std::array<Tvalue,N> temp;

  std::transform(values.begin(),values.end(),func.begin(),temp.begin(),[](auto a, auto b) { return a * b; } );

  Tvalue sum = std::accumulate(temp.begin(),temp.end(),zero,[](auto a, auto b) { return a + b; } );

  return sum;
}

template <class T, class Tvector> bool HexQuadratic3D<T,Tvector>::
  Bmatrix(const std::vector<Tvector> &totalCoordinates, const std::array<Tvector,NUM_NODES> &derivatives,  
  TMatrix<T,3,NUM_NODES> &B, T &jacobian)
{
// below #else is NOT correct
#if 1
  // coordinates of this element nodes
  TMatrix<T,NUM_NODES,3> coords = coordMatrix(totalCoordinates);

  // matrix for derivatives of shape functions
  TMatrix<T,NUM_NODES,3> func = shapeFunctionDerivativesMatrix(derivatives);

  TMatrix<T,3,NUM_NODES> funcT = func.transpose();

  // Jacobi matrix
  TMatrix<T,3,3> Jacobi = funcT * coords;

  // invert matrix
  TMatrix<T,3,3> J1 = +Jacobi;

  // check Jacobian value
  jacobian = Jacobi.lastDeterminant();

  if (std::abs(jacobian) < Jacobi.tolerance())
  {
    errorString_ = "Jacobian is close to zero when calculating derivatives on real coordinates";
    return false;
  }

  // fill matrix B[8 x 3] = J-1 * B^ (Connor & Brebbia, Rus edition, page 115, equation 9)
  // or J.J.Connor, C.A.Brebbia Finite Element Techniques for Fluid Flow, page 125, equation (i)
  B = J1 * funcT;

#else

  // coordinates of this element nodes
  TMatrix<T,NUM_NODES,3> coords = coordMatrix(totalCoordinates);

  // matrix for derivatives of shape functions
  TMatrix<T,NUM_NODES,3> func = shapeFunctionDerivativesMatrix(derivatives);

  // Jacobi matrix
  TMatrix<T,3,NUM_NODES> coordsT = coords.transpose();
  TMatrix<T,3,3> Jacobi = coordsT * func;
  // invert matrix
  TMatrix<T,3,3> J1 = +Jacobi;

  // check Jacobian value
  jacobian = Jacobi.lastDeterminant();
  if (std::abs(jacobian) < Jacobi.tolerance())
  {
    errorString_ = "Jacobian is close to zero when calculating derivatives on real coordinates";
    return false;
  }

  // fill matrix B[20 x 3] = J-1 * B^ (Connor & Brebbia, Rus edition, page 115, equation 9)
  // or J.J.Connor, C.A.Brebbia Finite Element Techniques for Fluid Flow, page 125, equation (i)
  TMatrix<T,3,NUM_NODES> funcT = func.transpose();
  B = J1 * funcT;

#endif

  return true;
}

template <class T, class Tvector> bool HexQuadratic3D<T,Tvector>::stiffnessMatrix(
  const std::vector<Tvector>& totalCoordinates,
  TMatrix<T,HexQuadratic3D<T,Tvector>::NUM_NODES,HexQuadratic3D<T,Tvector>::NUM_NODES>& matrix, const int power)
{
  // zero matrix
  matrix.fill(T0);

  // shape function derivatives
  for (int i = 0; i < GaussInt[power].numpoints; i++)
  {
    T e = static_cast<T>(GaussInt[power].knots[i]);
    T ew = static_cast<T>(GaussInt[power].weights[i]);
    for (int j = 0; j < GaussInt[power].numpoints; j++)
    {
      T n = static_cast<T>(GaussInt[power].knots[j]);
      T nw = static_cast<T>(GaussInt[power].weights[j]);
      for (int k = 0; k < GaussInt[power].numpoints; k++)
      {
        T c = static_cast<T>(GaussInt[power].knots[k]);
        T cw = static_cast<T>(GaussInt[power].weights[k]);

        // parms -1..+1
#if 1
        Tvector parms(e,n,c);
#else
        Tvector parms((e + T1) * T05,(n + T1) * T05,(c + T1) * T05);
#endif
        std::array<Tvector,NUM_NODES> funcderivatives;
        shapeDerivatives(parms,funcderivatives);

        TMatrix<T,3,NUM_NODES> B;
        T jacobian;
        bool ok = Bmatrix(totalCoordinates,funcderivatives,B,jacobian);
        if (!ok)
          return false;

        for (int l = 0; l < NUM_NODES; l++)
        {
          for (int m = 0; m < NUM_NODES; m++)
          {
            // add to integral
            T value = (B[0][l] * B[0][m] + B[1][l] * B[1][m] + B[2][l] * B[2][m]) * 
              std::abs(jacobian) * ew * nw * cw; 
            matrix[l][m] += value;
          }
        }
      }
    }
  }

  return true;
}

template <class T, class Tvector> bool HexQuadratic3D<T,Tvector>::massMatrix(
  const std::vector<Tvector>& totalCoordinates,
  TMatrix<T,HexQuadratic3D<T,Tvector>::NUM_NODES,HexQuadratic3D<T,Tvector>::NUM_NODES>& matrix, const int power)
{
  // zero matrix
  matrix.fill(T0);

  // shape function derivatives
  for (int i = 0; i < GaussInt[power].numpoints; i++)
  {
    T e = static_cast<T>(GaussInt[power].knots[i]);
    T ew = static_cast<T>(GaussInt[power].weights[i]);
    for (int j = 0; j < GaussInt[power].numpoints; j++)
    {
      T n = static_cast<T>(GaussInt[power].knots[j]);
      T nw = static_cast<T>(GaussInt[power].weights[j]);
      for (int k = 0; k < GaussInt[power].numpoints; k++)
      {
        T c = static_cast<T>(GaussInt[power].knots[k]);
        T cw = static_cast<T>(GaussInt[power].weights[k]);

        // parms -1..+1
#if 1
        Tvector parms(e,n,c);
#else
        Tvector parms((e + T1) * T05,(n + T1) * T05,(c + T1) * T05);
#endif
        std::array<T,NUM_NODES> func;
        shapeFunctions(parms,func);

        std::array<Tvector,NUM_NODES> funcderivatives;
        shapeDerivatives(parms,funcderivatives);

        TMatrix<T,3,NUM_NODES> B;
        T jacobian;
        bool ok = Bmatrix(totalCoordinates,funcderivatives,B,jacobian);
        if (!ok)
          return false;

        for (int l = 0; l < NUM_NODES; l++)
        {
          for (int m = 0; m < NUM_NODES; m++)
          {
            // add to integral
            T value = func[l] * func[m] * std::abs(jacobian) * ew * nw * cw;
            matrix[l][m] += value;
          }
        }
      }
    }
  }

  return true;
}

template <class T, class Tvector> T HexQuadratic3D<T,Tvector>::volume(
  const std::vector<Tvector>& totalCoordinates, const int power)
{
  T volume = T0;

  // shape function derivatives
  for (int i = 0; i < GaussInt[power].numpoints; i++)
  {
    T e = static_cast<T>(GaussInt[power].knots[i]);
    T ew = static_cast<T>(GaussInt[power].weights[i]);
    for (int j = 0; j < GaussInt[power].numpoints; j++)
    {
      T n = static_cast<T>(GaussInt[power].knots[j]);
      T nw = static_cast<T>(GaussInt[power].weights[j]);
      for (int k = 0; k < GaussInt[power].numpoints; k++)
      {
        T c = static_cast<T>(GaussInt[power].knots[k]);
        T cw = static_cast<T>(GaussInt[power].weights[k]);

        // parms -1..+1
#if 1
        Tvector parms(e,n,c);
#else
        Tvector parms((e + T1) * T05,(n + T1) * T05,(c + T1) * T05);
#endif
        std::array<T,NUM_NODES> func;
        shapeFunctions(parms,func);

        std::array<Tvector,NUM_NODES> funcderivatives;
        shapeDerivatives(parms,funcderivatives);

        TMatrix<T,3,NUM_NODES> B;
        T jacobian;
        bool ok = Bmatrix(totalCoordinates,funcderivatives,B,jacobian);
        if (!ok)
          return -1.0;

        volume += std::abs(jacobian) * ew * nw * cw;
      }
    }
  }

  return volume;
}

template <class T, class Tvector> void HexQuadratic3D<T,Tvector>::nodeMinMax(LINT& min, LINT& max)
{
  min = +std::numeric_limits<LINT>::max();
  max = -std::numeric_limits<LINT>::max();

  for (size_t j = 0; j < NUM_NODES; j++)
  {
    min = std::min(min,nodes_[j]);
    max = std::max(max,nodes_[j]);
  }
}

template <class T, class Tvector> Tvector HexQuadratic3D<T,Tvector>::coordinate(const Tvector& parms,
  const std::vector<Tvector>& totalCoordinates)
{
  // nodal coordinate values
  TMatrix<T,HexQuadratic3D<T,Tvector>::NUM_NODES,3> coords = coordMatrix(totalCoordinates);

  // shape functions
  std::array<T,NUM_NODES> func;
  shapeFunctions(parms,func);
  TMatrix<T,HexQuadratic3D<T,Tvector>::NUM_NODES,1> funcmatrix = shapeFunctionMatrix(func);

  // multiply
  TMatrix<T,3,1> res = coords.transpose() * funcmatrix;

  return Tvector(res[0][0],res[1][0],res[2][0]);
}

template <class T, class Tvector> T HexQuadratic3D<T,Tvector>::value(const Tvector& parms,
  const std::vector<T>& totalValues)
{
  // nodal values
  TMatrix<T,1,HexQuadratic3D<T,Tvector>::NUM_NODES> nvalues = values(totalValues);

  // shape functions
  std::array<T,NUM_NODES> func;
  shapeFunctions(parms,func);
  TMatrix<T,HexQuadratic3D<T,Tvector>::NUM_NODES,1> funcmatrix = shapeFunctionMatrix(func);

  // multiply
  TMatrix<T,1,1> res = nvalues * funcmatrix;

  return res[0][0];
}

template <class T, class Tvector> Tvector HexQuadratic3D<T,Tvector>::gradient(const Tvector& parms,
  const std::vector<Tvector> &totalCoordinates, const std::vector<T>& totalValues)
{
  // nodal values
  TMatrix<T,1,HexQuadratic3D<T,Tvector>::NUM_NODES> nvalues = values(totalValues);

  std::array<Tvector,NUM_NODES> funcderivatives;
  shapeDerivatives(parms,funcderivatives);

  TMatrix<T,3,NUM_NODES> B;
  T jacobian;
  bool ok = Bmatrix(totalCoordinates,funcderivatives,B,jacobian);

  TMatrix<T,3,1> d = B * nvalues.transpose();

  return Tvector(d[0][0],d[1][0],d[2][0]);
}
