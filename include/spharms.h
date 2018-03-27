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
  SphericalSurface(const natural_t rows, const natural_t cols);
  
  real_t get(const natural_t theta_n, const natural_t psi_m);
  void set(const natural_t theta_n, const natural_t psi_m, const real_t radius_nm);

  const natural_t rows;
  const natural_t cols;

private:
  std::vector<real_t> values_;
};

class SphericalHarmonics
{
  explicit SphericalHarmonics(const natural_t l_max);

  complex_t get(const natural_t l, const natural_t m);
  void set(const natural_t l, const natural_t m, const complex_t value);

  const natural_t l_max;
  
private:
  std::vector<complex_t> values_;
};

class Spharm
{
public:
  static SphericalHarmonics spharm_transform(const SphericalSurface& spherical_surface);
  static SphericalSurface ispharm_transform(const SphericalHarmonics& spherical_harmonics);
};

}
