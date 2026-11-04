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
     
#include <random>                         
#include "tbasics.h"
#include "tpoint.h"
#include "ttriangles.h"

using namespace std;
using namespace tcad;

/** Get random position on one of triangles. */
template <typename T> TPoint<T> getRandomTriPosition(TTriangles<T> &tris, std::mt19937& gen)
{
  // random tri number
  int N = int(tris.numFaces());
  int index = N * gen() / gen.max();
  LIMIT(index,0,N - 1);

  // barycentric coords
  T U = T(gen()) / T(gen.max());
  T V = T(gen()) / T(gen.max());
  T W = 1.0 - U - V;

  std::array<TPoint<T>,3> coords = tris.threeCorners(index);
  TPoint<T> p = coords[0] * U + coords[1] * V + coords[2] * W;

  return p;
}

/** Suggested by widely used Google search AI. */

// Avogadro number, # molecules in 1 mol
#define NA 6.022e23

// Boltzmann constant
#define KB 1.380649e-23 

/** "m" - mass of dry air molecule, default is for dry air, 0.02897 kg/mol. */
template <typename T> T m(T molar_mass = 0.02897)
{
  return molar_mass / NA;
}

/** Mean thermal molecular speed vn, must be something like 400 m/s. 
  temperature is in K. */
template <typename T> T meanThermalSpeed(T temperature, T molecular_mass)
{
  return sqrt(8.0 * KB * temperature / (PI10 * molecular_mass));
}

/** Max representative speed, mean thermal * 2, must be something like 400 m/s. 
  temperature is in K. */
template <typename T> T maxThermalSpeed(T temperature, T molecular_mass)
{
  return meanThermalSpeed(temperature,molecular_mass) * 3.0;
}

/** dt from fine mesh size. */
template <typename T> T dtFromCellSize(T cell_size, T temperature, T molecular_mass)
{
  return cell_size / maxThermalSpeed(temperature,molecular_mass);
}

/** dt from distance between rotor/stator. */
template <typename T> T dtFromBladeDistance(T d, T max_blade_speed)
{
  return 0.05 * d / max_blade_speed;
}

/** Max blade speed. R - rotor radius, RPM - # rotations per minute. */
template <typename T> T maxBladeSpeed(T R, T RPM)
{
  return PI20 * R / (60.0 / RPM);
}

/** Get number density n (# molecules per m^3). 1.0e-3 is target vacuum 
  pressure in Pa. */
template <typename T> T numberDensity(T temperature, T P = 1.0e-3)
{
  T n = P / (KB * temperature);
  return n;
}

/** Get total number of molecules inside simulation volume V. 
  1.0e-3 is target vacuum pressure in Pa. */
template <typename T> T Nphysical(T temperature, T volume, T P = 1.0e-3)
{
  T Np = numberDensity<T>(temperature,P) * volume;
  return Np;
}

/** Calculate specific weight Fnum (number of molecules in one particle). */
template <typename T> T Fnum(T temperature, T volume, int num_cells, T P = 1.0e-3, 
  int particles_per_cell = 10)
{
  T fnum = Nphysical<T>(temperature,volume,P) / T(num_cells * particles_per_cell);
  return fnum;
}

/** Nstep : how many particles to inject at any time step. 
  n - numberDensity(), 
  vn - meanThermalSpeed()
*/
template <typename T> int Nstep(T n, T inlet_area, T dt, T vn, T Fnum)
{
  int N = int(n * inlet_area * dt * vn / (4.0 * Fnum));
  LIMIT_MIN(N,1);
  return N;
}

/** Particle. */
template <typename T> struct Particle {
  TPoint<T> position;
  TPoint<T> velocity;
  bool active = false;
};

/** Generates a physically accurate random Maxwellian velocity vector for a new gas particle.
  tris are entry plane triangles
  gen must be initialised : mt19937 mt(time(nullptr)); 
  air molecular mass : 0.02896 kg/mol
  temperature : 20C
*/
template <typename T> void initialize_inlet_particle(Particle<T>& p, TTriangles<T> &tris,
  T v_bulk_z /* must be negative */, std::mt19937& gen,
  T temperature, T molecular_mass) 
{
  // Calculate the standard deviation (thermal spread) per axis
  T sigma = std::sqrt((KB * temperature) / molecular_mass);
    
//!!!!!!!
//T v = meanThermalSpeed(temperature,molecular_mass);

  // Set up the Gaussian (Normal) sampler
  std::normal_distribution<T> maxwell_axis(0.0,sigma);
    
  // Get random position on one of entry triangles
  p.position = getRandomTriPosition(tris,gen);

  // Sample continuous, independent random thermal velocities
  p.velocity.init(maxwell_axis(gen),maxwell_axis(gen),maxwell_axis(gen) + v_bulk_z); // Add directional suction flow if applicable

  // Start this particle
  p.active = true;
}

/** Executed immediately when your octree detects a particle-blade intersection. 
  blade_n is normalised. */
template <typename T> void apply_maxwellian_diffuse_reflection(Particle<T>& p, 
  TPoint<T> blade_v, TPoint<T> blade_n,
  T blade_temp, T molecular_mass, std::mt19937& gen) 
{
  T sigma = std::sqrt((KB * blade_temp) / molecular_mass);
    
  std::normal_distribution<T> gauss(0.0,sigma);
    
  // 1. Generate thermal velocities in a local surface coordinate system
  // Tangential components are pure Gaussians
  T v_tangent1 = gauss(gen);
  T v_tangent2 = gauss(gen);
    
  // Perpendicular component must be positively biased away from the wall.
  // A Maxwellian flux leaving a surface yields a Rayleigh distribution for the normal component:
  T r1 = gauss(gen);
  T r2 = gauss(gen);
  T v_normal = std::sqrt(r1*r1 + r2*r2); // Guaranteed positive, biased outward
    
  // 2. Map local surface velocities (v_tangent1, v_tangent2, v_normal) back to 3D Space (X, Y, Z)
  // using your octree's extracted triangle normal vector (nx, ny, nz)
  // (Assuming simple local alignment for this example; requires a coordinate basis rotation in production)
  TPoint<T> thermal_v;
  thermal_v.X = v_tangent1 + v_normal * blade_n.X;
  thermal_v.Y = v_tangent2 + v_normal * blade_n.Y;
  thermal_v.Z = gauss(gen) + v_normal * blade_n.Z; 

  // 3. Add the macroscopic moving velocity of the spinning blade (Momentum Transfer!)
  p.velocity = thermal_v + blade_v;
}
