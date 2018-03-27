/**
 * @file spharms.cpp
 * @author Sylvaus
 * @date Mon Jan 08 2018
 * @brief 
 *
 * Spherical Harmonics Helpers
 */

#include "spharms.h"

namespace hyperspharm
{


SphericalSurface::SphericalSurface(const natural_t rows, const natural_t cols) :
  rows(rows), cols(cols), values_(rows * cols)
{
}

real_t SphericalSurface::get(natural_t theta_n, const natural_t psi_m)
{
  return values_[rows * theta_n + psi_m];
}

void SphericalSurface::set(const natural_t theta_n, const natural_t psi_m, const real_t radius_nm)
{
  values_[rows * theta_n + psi_m] = radius_nm;
}



SphericalHarmonics::SphericalHarmonics(const natural_t l_max) :
  l_max(l_max), values_(((l_max + 2)*(l_max + 1))/2)
{
}

complex_t SphericalHarmonics::get(const natural_t l, const natural_t m)
{
  const size_t index = ((l+1) * l)/2 + m;
  if(index >= values_.size())
  {
    return {0, 0};
  }

  return values_[index];
}

void SphericalHarmonics::set(const natural_t l, const natural_t m, const complex_t value)
{
  const size_t index = ((l+1) * l)/2 + m;
  if(index >= values_.size())
  {
    return;
  }
  values_[index] = value;
}
}
