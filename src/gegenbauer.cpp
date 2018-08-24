/**
 * @file gegenbauer.cpp
 * @author Sylvaus
 * @date Tue July 08 2018
 * @brief
 *      Gegenbauer Polynomial helper functions
 *
 * Gegenbauer recurrence:
 *  C^l_m(x) = \frac{1}{l}[2x(l+m-1)C^{l-1}_m(x) - (l+2m-2)C^{l-2}_m(x)]
 */
#include "gegenbauer.h"

// TODO improve macros ((l + m) should be ((l) + (m))) or replace by inline if speed is same
#define GGALM(l, m) 2 * std::sqrt(static_cast<real_t>((l + m) * (l + m - 1)) \
                                  / static_cast<real_t>((l) * (l + (2 * m) - 1)))
#define GGBLM(l, m) - std::sqrt(static_cast<real_t>((l - 1) * (l + m) * (l + (2 * m) - 2)) \
                                 / static_cast<real_t>((l) * (l + m -2) * (l + (2 * m) - 1)))

#define GGA_WRITE_VALUE(values, l, m, value) values[m][l] = value
#define GGA_GET_VALUE(values, l, m) values[m][l]

namespace hyperspharm
{

std::vector<std::vector<GegenbauerPoly::coeff>> GegenbauerPoly::coeffs_ = {};

real_t GegenbauerPoly::get_normalized(const natural_t l, const natural_t m, const real_t x)
{
#ifndef NOCHECK
  if (std::abs(x) > 1.0)
  {
    throw std::invalid_argument( "Gegenbauer Polynomial: function domain is x in [-1, 1] " );
  }
  if (m < 1)
  {
    throw std::invalid_argument( "Gegenbauer Polynomial: m should be bigger than 0 " );
  }
#endif
  
  real_t N_0_m = N01;
  real_t N_1_m = 2.0 * x * N01;
  
  for (natural_t i = 2; i <= m; ++i)
  {
    N_0_m *= std::sqrt(1.0 + (1.0 / static_cast<real_t>(2 * i - 1)));
    N_1_m *= std::sqrt(1.0 + (3.0 / static_cast<real_t>(2 * i - 1)));
  }
  
  if (l == 0)
  {
    return N_0_m;
  }
    
  if (l == 1)
  {
    return N_1_m;
  }

  for (natural_t i = 2; i <= l; ++i)
  {
    const real_t tmp = N_1_m;

    N_1_m = (x * GGALM(i, m) * N_1_m) + (GGBLM(i, m) * N_0_m);
    N_0_m = tmp;
  }
  
  return N_1_m;
}

GegenbauerArray hyperspharm::GegenbauerPoly::get_norm_array(real_t norm_coeff, natural_t l_max, real_t x)
{
#ifndef NOCHECK
  if (std::abs(x) > 1.0)
  {
    throw std::invalid_argument( "Gegenbauer Polynomial: function domain is x in [-1, 1] " );
  }
  if (l_max < 1)
  {
    throw std::invalid_argument( "Gegenbauer Polynomial: l_max should be bigger than 0 " );
  }
#endif

  GegenbauerArray result(l_max);
  auto& values = result.values();
  compute_coefficients(l_max);
  
  real_t N_0_m = N01;
  real_t N_1_m = 2.0 * x * N01;

  
  for (natural_t m = 1; m <= l_max; ++m)
  {

    GGA_WRITE_VALUE(values, 0, m, N_0_m);
    GGA_WRITE_VALUE(values, 1, m, N_1_m);
    
    auto& coeffs_m = coeffs_[m];
    for (natural_t l = 2; l <= l_max; ++l)
    {
      GGA_WRITE_VALUE(values, l, m, 
                      (coeffs_m[l].a * x * GGA_GET_VALUE(values, l-1, m)) +
                      (coeffs_m[l].b * GGA_GET_VALUE(values, l-2, m)));
    }

    //Calculate the next N_0_m and N_1_m
    N_0_m *= std::sqrt(1.0 + (1.0 / static_cast<real_t>(2 * m + 1)));
    N_1_m *= std::sqrt(1.0 + (3.0 / static_cast<real_t>(2 * m + 1)));
  }
  
  for (auto& value : values)
  {
    for (auto& val : value)
    {
      val *= norm_coeff;
    }
  }
  
  return result;
}

void GegenbauerPoly::compute_coefficients(const natural_t l_max)
{
  if(coeffs_.size() > l_max) { return; }
  else { coeffs_.resize(l_max + 1); }

  for (natural_t m = 0; m <= l_max; ++m)
  {
    auto& coeffs_m = coeffs_[m];
    coeffs_m.resize(l_max + 1);
    for (natural_t l = 2; l <= l_max; ++l)
    {
      coeffs_m[l].a = GGALM(l, m);
      coeffs_m[l].b = GGBLM(l, m);
    }
  }
}



GegenbauerArray::GegenbauerArray(const natural_t l_max) :
  l_max_(l_max)
{
  values_.resize(l_max + 1);
  for (auto& value : values_)
  {
    value.resize(l_max + 1);
  }
}

real_t GegenbauerArray::get(const natural_t l, const natural_t m) const
{
  if(l > l_max_)
  {
    return 0.0;
  }

  return values_[m][l];
}

real_t GegenbauerArray::unsafe_get(const natural_t l, const natural_t m) const
{
  return values_[m][l];
}

void GegenbauerArray::set(const natural_t l, const natural_t m, const real_t x)
{
  if((m > l_max_) || (l > l_max_))
  {
    // TODO improve error message
    throw std::invalid_argument( "Gegenbauer array index out of range" );
  }
  
  values_[m][l] =  x;
}

void GegenbauerArray::unsafe_set(const natural_t l, const natural_t m, const real_t x)
{
  values_[m][l] =  x;
}

GegenbauerArray::GegenbauerArray(GegenbauerArray &&other) noexcept :
  l_max_(other.l_max_)
{
  values_ = std::move(other.values_);
}

GegenbauerArray &GegenbauerArray::operator=(GegenbauerArray &&other) noexcept
{
  if (this != &other)
  {
    l_max_ = other.l_max_;
    values_ = std::move(other.values_);
  }
  return *this;
}

natural_t GegenbauerArray::l_max() const 
{
  return l_max_;
}

std::vector<std::vector<real_t>>& GegenbauerArray::values()
{
  return values_;
}

}

