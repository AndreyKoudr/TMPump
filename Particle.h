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
#include "tmisc.h"
#include "ttriangles.h"

using namespace std;
using namespace tcad;

/** Get random position on one of triangles. */
template <typename T> TPoint<T> getRandomTriPosition(std::array<TPoint<T>,3> &coords, std::mt19937& gen)
{
  // barycentric coords
  T U = T(gen()) / T(gen.max());
  T V = T(gen()) / T(gen.max());
  T W = 1.0 - U - V;

  TPoint<T> p = coords[0] * U + coords[1] * V + coords[2] * W;

  return p;
}

/** Get random position on one of triangles. */
template <typename T> TPoint<T> getRandomTriPosition(TTriangles<T> &tris, std::mt19937& gen)
{
  // random tri number
  int N = int(tris.numFaces());
  int index = N * gen() / gen.max();
  LIMIT(index,0,N - 1);

  std::array<TPoint<T>,3> coords = tris.threeCorners(index);
  TPoint<T> p = getRandomTriPosition<T>(coords,gen);

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

/** Max representative speed, mean thermal * 3, must be something like 400 m/s. 
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
template <typename T> T bladeSpeed(T R, T RPM)
{
  return PI20 * R * RPM / 60.0;
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
  TPoint<T> oldposition;   // previous position
  TPoint<T> position;      // current position
  TPoint<T> velocity; // velocity
  bool active = false;
};

/** Generates a physically accurate random Maxwellian velocity vector for a new gas particle.
  gen must be initialised : mt19937 mt(time(nullptr)); 
  air molecular mass : 0.02896 kg/mol / Avogadro number
  temperature : Kelvins
  active is NOT SET.
*/
template <typename T> void setParticleVelocity(Particle<T>& p, 
  T v_bulk_z /* must be negative if not zero */, std::mt19937& gen,
  T temperature, T molecular_mass) 
{
  // Calculate the standard deviation (thermal spread) per axis
  T sigma = std::sqrt((KB * temperature) / molecular_mass);
    
  // Set up the Gaussian (Normal) sampler
  std::normal_distribution<T> maxwell_axis(0.0,sigma);
    
  // Sample continuous, independent random thermal velocities
  p.velocity.X = maxwell_axis(gen);
  p.velocity.Y = maxwell_axis(gen);
  p.velocity.Z = maxwell_axis(gen) + v_bulk_z;

  // there is something deeply wrong with this init()
  // p.velocity.init(maxwell_axis(gen),maxwell_axis(gen),maxwell_axis(gen) + v_bulk_z); // Add directional suction flow if applicable
}

#if 1

/** Executed immediately when your octree detects a particle-blade intersection. 
  blade_n is normalised. */
template <typename T> void applyParticleDiffuseReflection(Particle<T>& p, 
  TPoint<T> blade_v, TPoint<T> blade_n,
  T blade_temp, T molecular_mass, std::mt19937& gen) 
{
  // Most probable thermal speed: sqrt(2 * k_B * T / m)
  const T v_mp = std::sqrt((2.0 * KB * blade_temp) / molecular_mass);
  const T sigma_t = v_mp / std::sqrt(2.0); // Standard deviation for tangential components

  // Distributions (static or passed in to avoid re-initialization overhead)
  std::uniform_real_distribution<T> uniform_dist(0.0,1.0);
  std::normal_distribution<T> normal_dist(0.0,1.0);

  // 1. Sample local velocity components in the blade's frame
  // Shift (1.0 - U) so argument to std::log is strictly in (0.0, 1.0], preventing log(0) = -inf
  const T u1 = 1.0 - uniform_dist(gen); 
  const T v_n  = v_mp * std::sqrt(-std::log(u1));
  const T v_t1 = sigma_t * normal_dist(gen);
  const T v_t2 = sigma_t * normal_dist(gen);

  // 2. Transform local thermal velocity to global coordinates
  TPoint<T> t1,t2;
  buildOrthonormalBasis<T>(blade_n,t1,t2);
  const TPoint<T> thermal_v = (blade_n * v_n) + (t1 * v_t1) + (t2 * v_t2);

  // 3. Add the macroscopic moving velocity of the spinning blade (Momentum Transfer!)
  p.velocity = thermal_v + blade_v;

  // under no circumstances a particle can move inside body
  assert(thermal_v * blade_n >= 0.0);
  //if (thermal_v * blade_n < 0.0) <- incorrect, ity must be relative to blade
  //  p.velocity = -p.velocity;
}

#else

/** Executed immediately when your octree detects a particle-blade intersection. 
  blade_n is normalised. */
template <typename T> void applyParticleDiffuseReflection(Particle<T>& p, 
  TPoint<T> blade_v, TPoint<T> blade_n,
  T blade_temp, T molecular_mass, std::mt19937& gen) 
{
  T sigma = std::sqrt((KB * blade_temp) / molecular_mass);
  std::normal_distribution<T> gauss(0.0,sigma);

#if 1
  // 2. Generate a random diffuse cosine-distributed reflection vector relative to the normal
  T r1 = gauss(gen);
  T r2 = gauss(gen);
  T r3 = std::abs(gauss(gen)); // Always move away from the inward normal
  // in local coords
  TPoint<T> thermal_v(r1,r2,r3);
#else
    
  // 1. Generate thermal velocities in a local surface coordinate system
  // Tangential components are pure Gaussians
  T v_tangent1 = gauss(gen);
  T v_tangent2 = gauss(gen);
    
  // Perpendicular component must be positively biased away from the wall.
  // A Maxwellian flux leaving a surface yields a Rayleigh distribution for the normal component:
  T r1 = gauss(gen);
  T r2 = gauss(gen);
  T v_normal = std::sqrt(r1*r1 + r2*r2); // Guaranteed positive, biased outward
  // ... another way is to
  // double r3 (i.e. v_normal) = std::abs(gauss(gen)); // Always move away from the inward normal

  // there are three component in local coordinate system x,y,z : v_tangent1,v_tangent2 and v_normal 
    
  // 2. Map local surface velocities (v_tangent1, v_tangent2, v_normal) back to 3D Space (X, Y, Z)
  // using your octree's extracted triangle normal vector (nx, ny, nz)
  // (Assuming simple local alignment for this example; requires a coordinate basis rotation in production)
  // in local coords
  TPoint<T> thermal_v(v_tangent1,v_tangent2,v_normal);
#endif

#if 1
  // we need to rotate it into real normal coordinates
  TTransform<T> t;
  t.Rotate(+(TPoint<T>(0.0,0.0,1.0) ^ blade_n),blade_n < TPoint<T>(0.0,0.0,1.0));
  thermal_v = t.applyTransform(thermal_v);
#else
  TPoint<T> thermal_v;
  thermal_v.X = v_tangent1 + v_normal * blade_n.X;
  thermal_v.Y = v_tangent2 + v_normal * blade_n.Y;
  thermal_v.Z = gauss(gen) + v_normal * blade_n.Z; 
#endif

  // 3. Add the macroscopic moving velocity of the spinning blade (Momentum Transfer!)
  p.velocity = thermal_v + blade_v;

  // under no circumstances a particle can move inside body
  if (p.velocity * blade_n < 0.0)
    p.velocity = -p.velocity;
}

#endif


#if 0

// Diffuse reflection off a moving blade surface
void handle_blade_collision(Particle& p, double blade_vx, double blade_vy, double blade_vz, 
                             double nx, double ny, double nz, double thermal_speed, std::mt19937& gen) 
{
    // 1. Transform particle velocity to the blade's rest frame
    double v_rel_x = p.vx - blade_vx;
    double v_rel_y = p.vy - blade_vy;
    double v_rel_z = p.vz - blade_vz;

    // 2. Generate a random diffuse cosine-distributed reflection vector relative to the normal
    std::normal_distribution<double> gauss(0.0, thermal_speed);
    double r1 = gauss(gen);
    double r2 = gauss(gen);
    double r3 = std::abs(gauss(gen)); // Always move away from the inward normal

    // 3. Simple hemispherical scattering mapped to the surface normal (nx, ny, nz)
    // (For complex geometry, you would use a proper coordinate rotation matrix here)
    double v_scatter_x = r1 + r3 * nx;
    double v_scatter_y = r2 + r3 * ny;
    double v_scatter_z = r2 + r3 * nz;

    // 4. Transform back to the absolute laboratory frame
    p.vx = v_scatter_x + blade_vx;
    p.vy = v_scatter_y + blade_vy;
    p.vz = v_scatter_z + blade_vz;
}

#endif

/** Set random particle position on entry plane, tris are entry plane triangles,
  gen must be initialised : mt19937 mt(time(nullptr)); 
  active is NOT SET. */
template <typename T> void setParticleEntryPosition(Particle<T>& p, TTriangles<T> &tris,
  std::mt19937& gen)
{
  // Get random position on one of entry triangles
  p.position = getRandomTriPosition(tris,gen);
}

/** Find inactive particle; -1 if not found, its index if found. */
template <typename T> int findInactiveParticle(std::vector<Particle<T>> &particles)
{
  for (int i = 0; i < int(particles.size()); i++)
  {
    if (particles[i].active == false)
    {
      return i;
    }
  }
  return -1;
}

/** If inactive particle not found, particles are extended. */
template <typename T> int newParticleIndex(std::vector<Particle<T>> &particles, int increment = 1024)
{
  // find among inactive
  int index = findInactiveParticle(particles);
  if (index >= 0)
    return index;

  // extend
  index = int(particles.size());
  particles.resize(particles.size() + increment);

  return index;
}

/** Count active particles. */
template <typename T> int countParticles(std::vector<Particle<T>> &particles)
{
  int count = 0;

  for (auto &p : particles)
  {
    if (p.active)
      count++;
  }

  return count;
}

/** Rotate particle around Z by angle in degrees. */
template <typename T> void rotateParticle(Particle<T> &particle, T angledeg)
{
  TTransform<T> t;
  t.Rotate(TPoint<T>(0.0,0.0,1.0),angledeg * CPI);
  particle.position = t.applyTransform(particle.position);
  particle.oldposition = t.applyTransform(particle.oldposition);
  particle.velocity = t.applyTransform(particle.velocity);
}
