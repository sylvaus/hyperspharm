/**
 * @file spharms.h
 * @author Sylvaus
 * @date Mon Jan 08 2018
 * @brief 
 *
 * Spherical Harmonics Helpers
 */

#pragma once

#include <stdint.h>
#include <vector>
#include "types.h"

namespace hyperspharm
{

/**
 * @brief Sphere Container 
 * 
 * Contains NxM values corresponding to the radius_nm for the different angles theta_n (inclination), psi_m (azimuth)
 * 
 */
class SphericalSurface
{
public:
  const uint32_t N;
  const uint32_t M;
  
  SphericalSurface(const uint32_t N, const uint32_t M);
  
  std::vector<double>& operator[](const int index );
  
  void changeNorthPole(double value);
  void changeSouthPole(double value);

private:
  std::vector<std::vector<double>> values;
};

class SphericalHarmonics
{
  const uint32_t N;
  const uint32_t M;
  
  SphericalHarmonics(const uint32_t N, const uint32_t M);
  
  std::vector<double>& operator[](const int index );
private:
  std::vector<complex_t> values;
};

class Spharm
{
  static std::vector<double> an; 
  static std::vector<double> bnm;
  
  static SphericalHarmonics spharm_transform(const SphericalSurface& spherical_surface);
  static SphericalSurface ispharm_transform(const SphericalHarmonics& spherical_harmonics);
};

}
