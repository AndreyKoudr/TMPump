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

#include <algorithm>
#include <string>
#include "Vector.h"
#include "Plane.h"

// Matrix by multipler operation is always scalar double, not T
#define MULTIPLIER_SCALAR

/** Matrix with N rows and M columns */
template <typename T, size_t N, size_t M> class TMatrix : public std::array<std::array<T,M>,N>
{
public:

  // Constructor
  TMatrix()
  {
    fill(T(0.0));
  }

  // Constructor
  TMatrix(const std::array<T,M> &data)
  {
    (*this)[0] = data;
  }

  // Equal with tolerance
  bool operator == (const TMatrix<T,N,M> &other) const
  {
    for (int i = 0; i < N; i++)
    {
      for (int j = 0; j < M; j++)
      {
        if (std::abs((*this)[i][j] - other[i][j]) > tolerance_)
          return false;
      }
    }

    return true;
  }

  // Fill matrix with value
  void fill(const T value)
  {
    std::for_each(this->begin(),this->end(),
      [&](auto &v) { std::fill(v.begin(),v.end(),value); } );
  }

  // Transpose
  TMatrix<T,M,N> transpose() const
  {
    TMatrix<T,M,N> transposed;

    for (size_t i = 0; i < N; i++)
    {
      for (size_t j = 0; j < M; j++)
      {
        transposed[j][i] = (*this)[i][j];
      }
    }

    return transposed;
  }

  // Transpose in place (if N == M)
  void transposeSelf()
  {
    static_assert(N == M,"Self-transpose requires square matrix");

    TMatrix<T,M,N> transposed = transpose();

    *this = transposed;
  }

  // +
  TMatrix<T,N,M> operator + (const TMatrix<T,N,M> &other) const
  {
    TMatrix<T,N,M> result;

    for (int i = 0; i < N; i++)
    {
      for (int j = 0; j < M; j++)
      {
        result[i][j] = (*this)[i][j] + other[i][j];
      }
    }

    return result;
  }

  void operator += (const TMatrix<T,N,M> &other)
  {
    *this = *this + other;
  }

  // -
  TMatrix<T,N,M> operator - (const TMatrix<T,N,M> &other) const
  {
    TMatrix<T,N,M> result;

    for (int i = 0; i < N; i++)
    {
      for (int j = 0; j < M; j++)
      {
        result[i][j] = (*this)[i][j] - other[i][j];
      }
    }

    return result;
  }

  void operator -= (const TMatrix<T,N,M> &other)
  {
    *this = *this - other;
  }

  // * scalar
#ifdef MULTIPLIER_SCALAR

  TMatrix<T,N,M> operator * (const double scalar) const
  {
    TMatrix<T,N,M> result;

    for (int i = 0; i < N; i++)
    {
      for (int j = 0; j < M; j++)
      {
        result[i][j] = (*this)[i][j] * scalar;
      }
    }

    return result;
  }

  void operator *= (const double scalar)
  {
    *this = *this * scalar;
  }

#else

  TMatrix<T,N,M> operator * (const T mult) const
  {
    TMatrix<T,N,M> result;

    for (int i = 0; i < N; i++)
    {
      for (int j = 0; j < M; j++)
      {
        result[i][j] = (*this)[i][j] * mult;
      }
    }

    return result;
  }

  void operator *= (const T mult)
  {
    *this = *this * mult;
  }

#endif

  // * vector
  TVector<T> operator * (const TVector<T> point) const
  {
    static_assert(N == M && N == 4,"Transformation matrix must be square 4 x 4");

    TVector<T> result(0,0,0,0);

    for (int i = 0; i < N; i++)
    {
      for (int j = 0; j < M; j++)
      {
        result.XYZW[i] += (*this)[i][j] * point.XYZW[j];
      }
    }

    return result;
  }

  // determinant, available for 2 x 2, 3 x 3 and 4 x 4
  T determinant() const
  {
    static_assert(N == M,"Determinant is defined only for square matrices");

    T value = static_cast<T>(0.0);

    if (N == 2)
    {
      value = (*this)[0][0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0];
    } else if (N == 3)
    {
      value = 
        ((*this)[0][0]*((*this)[1][1]*(*this)[2][2]-(*this)[2][1]*(*this)[1][2])) - 
        ((*this)[0][1]*((*this)[1][0]*(*this)[2][2]-(*this)[1][2]*(*this)[2][0])) + 
        ((*this)[0][2]*((*this)[1][0]*(*this)[2][1]-(*this)[1][1]*(*this)[2][0]));
    } else if (N == 4)
    {
      value =
      (*this)[0][3]*(*this)[1][2]*(*this)[2][1]*(*this)[3][0] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2]*(*this)[3][0] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2]*(*this)[3][0]+
      (*this)[0][2]*(*this)[1][1]*(*this)[2][3]*(*this)[3][0] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0]*(*this)[3][1]+
      (*this)[0][3]*(*this)[1][0]*(*this)[2][2]*(*this)[3][1] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3]*(*this)[3][1] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3]*(*this)[3][1]+
      (*this)[0][3]*(*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1]*(*this)[3][2] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1]*(*this)[3][2]+
      (*this)[0][1]*(*this)[1][0]*(*this)[2][3]*(*this)[3][2] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0]*(*this)[3][3]+
      (*this)[0][2]*(*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1]*(*this)[3][3] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2]*(*this)[3][3] + (*this)[0][0]*(*this)[1][1]*(*this)[2][2]*(*this)[3][3];
    } else
    {
      assert(false && "Determinant is calculated for square matrices 2 x 2, 3 x 3 and 4 x 4 only");
    }

    return value;
  }

  // Inverse matrix, 2 x 2, 3 x 3 and 4 x 4 only
  TMatrix<T,N,M> operator + ()
  {
    static_assert(N == M,"Only square matrix can be inverted");

    determinant_ = determinant();
    assert(std::abs(determinant_) >= std::numeric_limits<T>::epsilon() * static_cast<T>(10.0));

    TMatrix<T,N,M> res;

    // failure
    if (std::abs(determinant_) < tolerance_)
    {
      errorString_ = "TMatrix inversion failed due to close-to-zero determinant";
      res.fill(T(0.0));
      return res;
    }

    T invdet = static_cast<T>(1.0) / determinant_;

    if (N == 2)
    {
      res[0][0] =  +(*this)[1][1];
      res[0][1] =  -(*this)[0][1];
      res[1][0] =  -(*this)[1][0];
      res[1][1] =  +(*this)[0][0];
    } else if (N == 3)
    {
      res[0][0] =  ((*this)[1][1]*(*this)[2][2]-(*this)[1][2]*(*this)[2][1]);
      res[0][1] = -((*this)[0][1]*(*this)[2][2]-(*this)[2][1]*(*this)[0][2]);
      res[0][2] =  ((*this)[0][1]*(*this)[1][2]-(*this)[0][2]*(*this)[1][1]);
      res[1][0] = -((*this)[1][0]*(*this)[2][2]-(*this)[2][0]*(*this)[1][2]);
      res[1][1] =  ((*this)[0][0]*(*this)[2][2]-(*this)[2][0]*(*this)[0][2]);
      res[1][2] = -((*this)[0][0]*(*this)[1][2]-(*this)[0][2]*(*this)[1][0]);
      res[2][0] =  ((*this)[1][0]*(*this)[2][1]-(*this)[2][0]*(*this)[1][1]);
      res[2][1] = -((*this)[0][0]*(*this)[2][1]-(*this)[0][1]*(*this)[2][0]);
      res[2][2] =  ((*this)[0][0]*(*this)[1][1]-(*this)[0][1]*(*this)[1][0]);
    } else if (N == 4)
    {
      res[0][0] = (*this)[1][2]*(*this)[2][3]*(*this)[3][1] - (*this)[1][3]*(*this)[2][2]*(*this)[3][1] + (*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[1][1]*(*this)[2][2]*(*this)[3][3];
      res[0][1] = (*this)[0][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][2]*(*this)[2][3]*(*this)[3][1] - (*this)[0][3]*(*this)[2][1]*(*this)[3][2] + (*this)[0][1]*(*this)[2][3]*(*this)[3][2] + (*this)[0][2]*(*this)[2][1]*(*this)[3][3] - (*this)[0][1]*(*this)[2][2]*(*this)[3][3];
      res[0][2] = (*this)[0][2]*(*this)[1][3]*(*this)[3][1] - (*this)[0][3]*(*this)[1][2]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[3][2] - (*this)[0][1]*(*this)[1][3]*(*this)[3][2] - (*this)[0][2]*(*this)[1][1]*(*this)[3][3] + (*this)[0][1]*(*this)[1][2]*(*this)[3][3];
      res[0][3] = (*this)[0][3]*(*this)[1][2]*(*this)[2][1] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2] + (*this)[0][2]*(*this)[1][1]*(*this)[2][3] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3];
      res[1][0] = (*this)[1][3]*(*this)[2][2]*(*this)[3][0] - (*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[1][3]*(*this)[2][0]*(*this)[3][2] + (*this)[1][0]*(*this)[2][3]*(*this)[3][2] + (*this)[1][2]*(*this)[2][0]*(*this)[3][3] - (*this)[1][0]*(*this)[2][2]*(*this)[3][3];
      res[1][1] = (*this)[0][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][3]*(*this)[2][2]*(*this)[3][0] + (*this)[0][3]*(*this)[2][0]*(*this)[3][2] - (*this)[0][0]*(*this)[2][3]*(*this)[3][2] - (*this)[0][2]*(*this)[2][0]*(*this)[3][3] + (*this)[0][0]*(*this)[2][2]*(*this)[3][3];
      res[1][2] = (*this)[0][3]*(*this)[1][2]*(*this)[3][0] - (*this)[0][2]*(*this)[1][3]*(*this)[3][0] - (*this)[0][3]*(*this)[1][0]*(*this)[3][2] + (*this)[0][0]*(*this)[1][3]*(*this)[3][2] + (*this)[0][2]*(*this)[1][0]*(*this)[3][3] - (*this)[0][0]*(*this)[1][2]*(*this)[3][3];
      res[1][3] = (*this)[0][2]*(*this)[1][3]*(*this)[2][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3];
      res[2][0] = (*this)[1][1]*(*this)[2][3]*(*this)[3][0] - (*this)[1][3]*(*this)[2][1]*(*this)[3][0] + (*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[1][0]*(*this)[2][1]*(*this)[3][3];
      res[2][1] = (*this)[0][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][1]*(*this)[2][3]*(*this)[3][0] - (*this)[0][3]*(*this)[2][0]*(*this)[3][1] + (*this)[0][0]*(*this)[2][3]*(*this)[3][1] + (*this)[0][1]*(*this)[2][0]*(*this)[3][3] - (*this)[0][0]*(*this)[2][1]*(*this)[3][3];
      res[2][2] = (*this)[0][1]*(*this)[1][3]*(*this)[3][0] - (*this)[0][3]*(*this)[1][1]*(*this)[3][0] + (*this)[0][3]*(*this)[1][0]*(*this)[3][1] - (*this)[0][0]*(*this)[1][3]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[3][3] + (*this)[0][0]*(*this)[1][1]*(*this)[3][3];
      res[2][3] = (*this)[0][3]*(*this)[1][1]*(*this)[2][0] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1] + (*this)[0][1]*(*this)[1][0]*(*this)[2][3] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3];
      res[3][0] = (*this)[1][2]*(*this)[2][1]*(*this)[3][0] - (*this)[1][1]*(*this)[2][2]*(*this)[3][0] - (*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[1][0]*(*this)[2][2]*(*this)[3][1] + (*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[1][0]*(*this)[2][1]*(*this)[3][2];
      res[3][1] = (*this)[0][1]*(*this)[2][2]*(*this)[3][0] - (*this)[0][2]*(*this)[2][1]*(*this)[3][0] + (*this)[0][2]*(*this)[2][0]*(*this)[3][1] - (*this)[0][0]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[2][0]*(*this)[3][2] + (*this)[0][0]*(*this)[2][1]*(*this)[3][2];
      res[3][2] = (*this)[0][2]*(*this)[1][1]*(*this)[3][0] - (*this)[0][1]*(*this)[1][2]*(*this)[3][0] - (*this)[0][2]*(*this)[1][0]*(*this)[3][1] + (*this)[0][0]*(*this)[1][2]*(*this)[3][1] + (*this)[0][1]*(*this)[1][0]*(*this)[3][2] - (*this)[0][0]*(*this)[1][1]*(*this)[3][2];
      res[3][3] = (*this)[0][1]*(*this)[1][2]*(*this)[2][0] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2] + (*this)[0][0]*(*this)[1][1]*(*this)[2][2];
    } else
    {
      assert(false && "Only square matrices 2 x 2, 3 x 3 and 4 x 4 can be inverted");
      res.fill(static_cast<T>(0.0));
    }

    res *= invdet;

    return res;
  }

  // Is it identity matrix?
  bool isIdentity() const
  { 
    if (N != M)
      return false;

    for (int i = 0; i < N; i++)
    {
      for (int j = 0; j < M; j++)
      {
        if (i == j)
        {
          if (std::abs((*this)[i][j] - static_cast<T>(1.0)) > tolerance_)
            return false;
        } else
        {
          if (std::abs((*this)[i][j]) > tolerance_)
            return false;
        }
      }
    }

    return true;
  }

  // Set identity matrix (only for square matrix)
  void setIdentity()
  {
    static_assert(N == M,"Setting identity requires square matrix");

    fill(static_cast<T>(0.0));

    setDiag(static_cast<T>(1.0));
  }

  // Set equal diagonal elements
  void setDiag(const T value)
  {
    static_assert(N == M,"Setting diagonal requires square matrix");

    for (int i = 0; i < N; i++)
    {
      for (int j = 0; j < M; j++)
      {
        if (i == j)
          (*this)[i][j] = value;
      }
    }
  }

  //===== Transformations =====
  
   // Set translation matrix, 4 x 4 only
  void setTranslation(T x, T y, T z)
  {
    static_assert(N == M && N == 4,"Transformation matrix must be square 4 x 4");

    fill(static_cast<T>(0.0));

    (*this)[0][0] = static_cast<T>(1.0);
    (*this)[1][1] = static_cast<T>(1.0);
    (*this)[2][2] = static_cast<T>(1.0);
    (*this)[3][3] = static_cast<T>(1.0);

    (*this)[0][3] = x;
    (*this)[1][3] = y;
    (*this)[2][3] = z;
  }

  // Set scaling matrix, 4 x 4 only
  void setScaling(T coefx, T coefy, T coefz)
  {
    static_assert(N == M && N == 4,"Transformation matrix must be square 4 x 4");

    fill(static_cast<T>(0.0));

    (*this)[0][0] = coefx;
    (*this)[1][1] = coefy;
    (*this)[2][2] = coefz;
    (*this)[3][3] = static_cast<T>(1.0);
  }

  // Set rotation, angles.X - different, angles.Y - course, angles.Z - roll, 4 x 4 only
  void setRotation(TVector<T> angles)
  {
    static_assert(N == M && N == 4,"Transformation matrix must be square 4 x 4");

    fill(static_cast<T>(0.0));

	  T sincourse = sin(-angles.X);
	  T coscourse = cos(-angles.X);

	  T sindifferent = sin(-angles.Y);
	  T cosdifferent = cos(-angles.Y);
  
	  T sinroll = sin(-angles.Z);
	  T cosroll = cos(-angles.Z);

	  T SBCA = sindifferent * coscourse;
	  T SBSA = sindifferent * sincourse;

    (*this)[0][0] = cosdifferent * coscourse;
    (*this)[0][1] = sindifferent;
    (*this)[0][2] = cosdifferent * sincourse;

    (*this)[1][0] = -SBCA * cosroll + sincourse * sinroll;
    (*this)[1][1] = cosdifferent * cosroll;
    (*this)[1][2] = -SBSA * cosroll - sinroll * coscourse;

    (*this)[2][0] = -SBCA * sinroll - sincourse * cosroll;
    (*this)[2][1] = cosdifferent * sinroll;
    (*this)[2][2] = -SBSA * sinroll+coscourse * cosroll;

    (*this)[3][3] = static_cast<T>(1.0);
  }

  // Set rotation matrix; WORKS WRONG for theta > 180, use -(360 - theta) !!!, angle in radians, 4 x 4 only
  void setRotation(T Ux, T Uy, T Uz, T theta)
  {
    static_assert(N == M && N == 4,"Transformation matrix must be square 4 x 4");

    fill(static_cast<T>(0.0));

    T c = cos(theta);
    T s = sin(theta);
    T c1 = 1.0 - c;

    (*this)[0][0] = Ux * Ux + c * (1.0 - Ux * Ux);
    (*this)[1][1] = Uy * Uy + c * (1.0 - Uy * Uy);
    (*this)[2][2] = Uz * Uz + c * (1.0 - Uz * Uz);
    (*this)[3][3] = 1.0;
    
    (*this)[0][1] = Ux * Uy * c1 - Uz * s;
    (*this)[0][2] = Uz * Ux * c1 + Uy * s;
    (*this)[1][0] = Ux * Uy * c1 + Uz * s;
    (*this)[1][2] = Uy * Uz * c1 - Ux * s;
    (*this)[2][0] = Uz * Ux * c1 - Uy * s;
    (*this)[2][1] = Uy * Uz * c1 + Ux * s;
  }

  // Set scaling + translation matrix, 4 x 4 only
  void setScalingTranslation(T sx, T sy, T sz, T x, T y, T z)
  {
    static_assert(N == M && N == 4,"Transformation matrix must be square 4 x 4");

    fill(static_cast<T>(0.0));

    (*this)[0][0] = sx;
    (*this)[1][1] = sy;
    (*this)[2][2] = sz;
    (*this)[3][3] = 1.0;

    (*this)[0][3] = x;
    (*this)[1][3] = y;
    (*this)[2][3] = z;
  }

  // Set perspective projection matrix, X points forward, Y - up, Z - to the right, 4 x 4 only;
  // d - distance to screen
  void setPerspective(T d)
  {
    static_assert(N == M && N == 4,"Transformation matrix must be square 4 x 4");

    fill(static_cast<T>(0.0));

    (*this)[0][0] = 1.0;
    (*this)[1][1] = 1.0;
    (*this)[2][2] = 1.0;
    (*this)[3][0] = 1.0 / d;
  }

  // Cast a shadow on plane from direction, 4 x 4 only
  void setShadow(Plane<T> *plane, TVector<T> *dir)
  {
    (*this)[0][0] = plane->Normal.Y * dir->Y + plane->Normal.Z * dir->Z;
    (*this)[0][1] = -(plane->Normal.X * dir->Y);
    (*this)[0][2] = -(plane->Normal.X * dir->Z);
    (*this)[0][3] = 0.0;
    (*this)[1][0] = -(plane->Normal.Y * dir->X);
    (*this)[1][1] = plane->Normal.X * dir->X + plane->Normal.Z * dir->Z;
    (*this)[1][2] = -(plane->Normal.Y * dir->Z);
    (*this)[1][3] = 0.0;
    (*this)[2][0] = -(plane->Normal.Z * dir->X);
    (*this)[2][1] = -(plane->Normal.Z * dir->Y);
    (*this)[2][2] = plane->Normal.X * dir->X + plane->Normal.Y * dir->Y;
    (*this)[2][3] = 0.0;
    (*this)[3][0] = -(plane->Constant * dir->X);
    (*this)[3][1] = -(plane->Constant * dir->Y);
    (*this)[3][2] = -(plane->Constant * dir->Z);
    (*this)[3][3] = plane->Normal.X * dir->X + plane->Normal.Y * dir->Y + plane->Normal.Z * dir->Z;
  }

  //===== Tensor stuff =====

  // Frobenius inner product ":" (result is a single T)
  T innerProduct(const TMatrix<T,N,M> &other) const
  {
    T result = T(0.0);

    for (int i = 0; i < N; i++)
    {
      for (int j = 0; j < M; j++)
      {
        result += (*this)[i][j] * other[i][j];
      }
    }

    return result;
  }

  //                T
  // sym = 0.5(M + M )
  TMatrix<T,N,M> sym() const
  {
    static_assert(N == M,"Tensor sym requires a square matrix");

    TMatrix<T,N,M> result = (*this + transpose()) * 0.5;

    return result;
  }

  // Trace (sum of diagonal elements)
  T trace() const
  {
    static_assert(N == M,"Tensor trace requires a square matrix");

    T result = T(0.0);

    for (int i = 0; i < N; i++)
    {
      result += (*this)[i][i];
    }

    return result;
  }

  //===== Tolerance =====

  // Tolerance
  static T tolerance()
  {
    return tolerance_;
  }

  // Set tolerance
  static void setTolerance(const T tolerance)
  {
    tolerance_ = tolerance;
  }

  // Set default tolerance
  static void setDefaultTolerance()
  {
    tolerance_ = std::numeric_limits<T>::epsilon() * static_cast<T>(100.0);
  }

  // Last determinant calculated during matrix inversion
  T lastDeterminant() const
  {
    return determinant_;
  }

  // Return last error string
  static std::string errorString()
  {
    return errorString_;
  }

  // Debug output
  void DebugMatrixOutput(char *name, int numdigits)
  {
    debugString(std::string(name) + std::string(" ") + std::string(int(N)) + std::string(" x ") + std::string(int(M)));
    for (int i = 0; i < N; i++)
    {
      std::string s;
      for (int j = 0; j < M; j++)
      {
        s = s + std::string(T((*this)[i][j]),numdigits) + std::string(" ");
      }
      ::debugString(s);
   }
  }

private:
  // determinant after last inversion
  T determinant_;
//!!!!!!!!!  T determinant_ = T(1.0);

  // tolerance
  static T tolerance_;

  // last error string
  static std::string errorString_;
};

template <typename T, size_t N, size_t M>
  T TMatrix<T,N,M>::tolerance_ = std::numeric_limits<T>::epsilon() * static_cast<T>(100.0);

template <typename T, size_t N, size_t M>
  std::string TMatrix<T,N,M>::errorString_;

template <typename T, size_t I, size_t K, size_t J> TMatrix<T,I,J> operator * 
  (const TMatrix<T,I,K> &m0, const TMatrix<T,K,J> &m1)
{
  TMatrix<T,I,J> result;

  for (int i = 0; i < I; i++)
  {
    for (int j = 0; j < J; j++)
    {
      result[i][j] = static_cast<T>(0.0);
      for (int k = 0; k < K; k++)
      { 
        result[i][j] += m0[i][k] * m1[k][j];
      }
    }
  }

  return result;
}
