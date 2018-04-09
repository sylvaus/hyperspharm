/**
 * @file spharms.h
 * @author Sylvaus
 * @date Mon Jan 08 2018
 * @brief 
 *
 * Spherical Harmonics Helpers
 */

#pragma once

#include <cstdint>
#include <algorithm>
#include <vector>
#include <iomanip>
#include "fft.h"
#include "utils.h"
#include "types.h"
#include "legendre.h"

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
  
  real_t get(const natural_t theta_n, const natural_t psi_m) const;
  void set(const natural_t theta_n, const natural_t psi_m, const real_t radius_nm);

  natural_t rows() const;
  natural_t cols() const;

  std::vector<complex_t> get_psi_array(const natural_t theta_n) const;

  std::string to_string();
private:
  natural_t rows_;
  natural_t cols_;
  std::vector<real_t> values_;
};

class SphericalHarmonics
{
public:
  explicit SphericalHarmonics(const natural_t l_max);

  complex_t get(const natural_t l, const natural_t m) const;
  void set(const natural_t l, const natural_t m, const complex_t value);

  natural_t l_max() const;

  std::string to_string();
private:
  natural_t l_max_;
  std::vector<complex_t> values_;
};

class Spharm
{
public:
  static SphericalHarmonics spharm_transform(const SphericalSurface& spherical_surface);
  static SphericalSurface ispharm_transform(const SphericalHarmonics& spherical_harmonics);
};

}
