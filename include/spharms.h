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
#include "utils.h"
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
  SphericalSurface(const natural_t nrows, const natural_t ncols);
  
  std::vector<real_t>& operator[](const int index );
  
  natural_t nrows();
  natural_t ncols();
  
  void changeNorthPole(const real_t value);
  void changeSouthPole(const real_t value);

private:
  natural_t nrows_;
  natural_t ncols_;
  std::vector<std::vector<real_t>> values_;
};

class SphericalHarmonics
{
  SphericalHarmonics(const natural_t nmax, const natural_t mmax);
  
  std::vector<real_t>& operator[](const natural_t index );
  
  natural_t nmax();
  natural_t mmax();
  
private:
  natural_t nmax_;
  natural_t mmax_;
  std::vector<complex_t> values_;
};

class Spharm
{
public:
  static void compute_clmk(const natural_t l_max);
  
  
  
  static SphericalHarmonics spharm_transform(const SphericalSurface& spherical_surface);
  static SphericalSurface ispharm_transform(const SphericalHarmonics& spherical_harmonics);
};

}
