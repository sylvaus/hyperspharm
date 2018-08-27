/**
 * @file legendre.cpp
 * @author Sylvaus
 * @date Mon Jan 08 2018
 * @brief 
 *
 * Legendre Polynomial Helpers
 * Based on following paper for Associated Legendre Polynomial: http://www.scielo.org.co/pdf/racefn/v37n145/v37n145a09.pdf
 * Based on following page for fully normalized Legendre Polynomial: http://mitgcm.org/~mlosch/geoidcookbook/node11.html
 */

#include "legendre.h"

namespace hyperspharm
{

const real_t LegendrePoly::SQRT_4_PI = sqrt(4.0 * M_PI);
const real_t LegendrePoly::SPHARM_NORM = 1.0 / SQRT_4_PI;
const real_t LegendrePoly::FULLY_NORM = 1.0 / std::sqrt(2.0);

std::vector<std::vector<LegendrePoly::coeff>> LegendrePoly::coeffs_ = {};

real_t LegendrePoly::get(const natural_t l, const real_t x)
{
  return get_associated(l, 0, x);
}

real_t LegendrePoly::get_associated(const natural_t l, 
                                    const integer_t m, 
                                    const real_t x)
{  
#ifndef NOCHECK
  if (static_cast<natural_t>(std::abs(m)) > l)
  {
    throw std::invalid_argument( "Associated Legendre Polynomial: abs(m) must be smaller or equal to l" );
  }

  if (std::abs(x) > 1.0)
  {
    throw std::invalid_argument( "Associated Legendre Polynomial: function domain is x in [-1, 1] " );
  }
#else
  const auto abs_m = static_cast<natural_t >(std::abs(m));
#endif

  if (almost_equal(std::abs(x), static_cast<real_t>(1.0), static_cast<real_t>(0.000000001)))
  {
    if (m != 0) {return 0;}
    if (x > 0.0) {return 1.0;}
    else {return minus_one_power(l);}
  }
 
  const real_t sqrt_1_x2 = sqrt(static_cast<real_t>(1.0) - x*x);
  const real_t half_sqrt_1_x2 = sqrt_1_x2 * 0.5;
  const real_t inv_half_sqrt_1_x2 = x / half_sqrt_1_x2;
  
  if (l == 0)
  {
    return 1;
  }
  else if (l == 1)
  {
    switch (m)
    {
      case -1:
        return 0.5 * sqrt_1_x2;
      case 0:
        return x;
      case 1:
        return -sqrt_1_x2;
      default:
      {/*Do nothing*/}
    }
  }

  real_t p_l__m = minus_one_power(l) * std::pow(half_sqrt_1_x2, l) / Factorial::get(l); // P_l^{-l}
  if (m == -static_cast<integer_t>(l)) {return p_l__m;}
  real_t p_l_1_m = -p_l__m * static_cast<real_t>(l) * inv_half_sqrt_1_x2; // P_l^{1-l}

  const auto abs_m = std::abs(m);
  for (integer_t k = 2 - l; k <= -abs_m; k++)
  {
    real_t tmp = p_l_1_m;
    p_l_1_m = (static_cast<real_t>(k-1) * (inv_half_sqrt_1_x2) * p_l_1_m)
              - (static_cast<real_t>(l - (k-2)) * static_cast<real_t>(l + (k - 1)) * p_l__m);
    p_l__m = tmp;
  }

  if (m == abs_m)
  {
    p_l_1_m = minus_one_power(abs_m) * (Factorial::get(l + abs_m)/Factorial::get(l - abs_m)) * p_l_1_m;
  }

  return is_even(m) ? p_l_1_m : -p_l_1_m;
}


real_t LegendrePoly::get_fully_normalized(const natural_t l, 
                                          const integer_t m, 
                                          const real_t x)
{
#ifndef NOCHECK
  const auto abs_m = static_cast<natural_t>(std::abs(m));
  if (abs_m > l)
  {
    throw std::invalid_argument( "Normalized Legendre Polynomial: abs(m) must be smaller or equal to l" );
  }

  if (std::abs(x) > 1.0)
  {
    throw std::invalid_argument( "Normalized Legendre Polynomial: function domain is x in [-1, 1] " );
  }
#else
  const auto abs_m = static_cast<natural_t >(std::abs(m));
#endif

  // Compute N_m^m
  real_t n_m_m = 1;
  if (m != 0)
  {
    const real_t poly = (1.0 - (x * x));
    for (natural_t i = 1; i <= abs_m; i++)
    {
      n_m_m *= std::sqrt(poly *
                         (2.0 * static_cast<real_t>(i) + 1.0) /
                         (2.0 * static_cast<real_t>(i)));
    }
    n_m_m = minus_one_power(abs_m) * (n_m_m);
  }
  
  if (l == abs_m) {return n_m_m;}
  if (l == -abs_m) 
  {
    if (is_even(abs_m)) {return n_m_m;}
    else {return -n_m_m;}
  }
  
  // Compute N_{m + 1}^m
  const auto abs_m_real = static_cast<real_t>(abs_m);
  
  real_t n_m_1_m = x * std::sqrt((2.0 * abs_m_real) + 3.0) * n_m_m;
  for (natural_t i = (abs_m + 2); i <= l; i++)
  {
    const real_t coeff_n_m_1_m = std::sqrt(static_cast<real_t>((4 * i * i) - 1) / // (2*i - 1) * (2*i + 1) = (2*i)²  - 1 = (4 * i²) -1
                                           static_cast<real_t>((i - abs_m) * (i + abs_m))
    );
    const real_t coeff_n_m_m = std::sqrt(static_cast<real_t>(((2 * i) + 1) * (i - abs_m - 1) * (i + abs_m - 1))/
                                         static_cast<real_t>((i - abs_m) * (i + abs_m) * ((2 * i) - 3))
    );
    const real_t tmp = n_m_1_m;
    n_m_1_m = (coeff_n_m_1_m * x * n_m_1_m) - 
              (coeff_n_m_m * n_m_m);
    n_m_m = tmp;
  }
  
  if (m < 0) 
  {
    if (is_even(abs_m)) {return n_m_1_m;}
    else {return -n_m_1_m;}
  }
  else
  {
    return n_m_1_m;
  }
}

real_t LegendrePoly::get_spharm_normalized(const natural_t l, 
                                          const integer_t m, 
                                          const real_t x)
{
  return get_fully_normalized(l, m, x) / SQRT_4_PI;
}


NormalizedLegendreArray LegendrePoly::get_norm_array(const real_t normalization_coeff,
                                                     const natural_t l_max,
                                                     const real_t x)
{
  NormalizedLegendreArray result(l_max);

  // Compute N_m^m
  result.values_[0][0] = 1.0;
  real_t pi_i = 1.0;
  if (l_max != 0)
  {
    const real_t poly_sqrt = -std::sqrt(static_cast<real_t>(1.0) - (x * x));
    for (natural_t i = 1; i <= l_max; i++)
    {
      pi_i = poly_sqrt * (pi_i * std::sqrt(static_cast<real_t>(2*i + 1) / static_cast<real_t>(2*i)));
      result.values_[i][i] = pi_i;
    }
  }

  compute_coefficients(l_max);
  for (natural_t m = 0; m < l_max; ++m)
  {
    auto& values_m = result.values_[m];
    values_m[m + 1] = x * std::sqrt(static_cast<real_t>(m * 2 + 3)) * values_m[m];

    auto& coeffs_m = coeffs_[m];
    for (natural_t l = (m + 2); l <= l_max; ++l)
    {
      values_m[l] = (coeffs_m[l].alm * x * values_m[l-1]) -
                    (coeffs_m[l].blm * values_m[l-2]);
    }
  }

  for (auto& values : result.values_)
  {
    for (auto& value : values)
    {
      value *= normalization_coeff;
    }
  }

  return result;
}

NormalizedLegendreArray LegendrePoly::get_fully_norm_array(const natural_t l_max, const real_t x)
{
  return get_norm_array(FULLY_NORM, l_max, x);
}

NormalizedLegendreArray LegendrePoly::get_sph_norm_array(const natural_t l_max, const real_t x)
{
  return get_norm_array(SPHARM_NORM, l_max, x);
}

void LegendrePoly::compute_coefficients(const natural_t l_max)
{
  if(coeffs_.size() > l_max) { return; }
  else { coeffs_.resize(l_max + 1); }

  for (natural_t m = 0; m <= l_max; ++m)
  {
    auto& coeffs_m = coeffs_[m];
    coeffs_m.resize(l_max + 1);
    for (natural_t l = (m + 2); l <= l_max; ++l)
    {
      coeffs_m[l].alm = std::sqrt(
          static_cast<real_t>((4 * l * l) - 1) / // (2*i - 1) * (2*i + 1) = (2*i)²  - 1 = (4 * i²) -1
          static_cast<real_t>((l - m) * (m + l))
      );
      coeffs_m[l].blm = std::sqrt(static_cast<real_t>(((2 * l) + 1) * (l - m - 1) * (l + m - 1)) /
                                  static_cast<real_t>((l - m) * (l + m) * ((2 * l) - 3))
      );

    }
  }
}

NormalizedLegendreArray::NormalizedLegendreArray(const natural_t l_max) :
  l_max_(l_max)
{
  values_.resize(l_max + 1);
  for (auto& value : values_)
  {
    value.resize(l_max + 1);
  }
}

real_t NormalizedLegendreArray::get(const natural_t l, const integer_t m) const
{
  if((static_cast<natural_t>(std::abs(m)) > l )|| l > l_max_)
  {
    return 0;
  }

  return ((m >= 0) || is_even(m)) ? values_[m][l] : -values_[m][l];
}

real_t NormalizedLegendreArray::unsafe_get(const natural_t l, const natural_t m) const
{
  return values_[m][l];
}

void NormalizedLegendreArray::set(const natural_t l, const integer_t m, const real_t x)
{
  if((static_cast<natural_t>(std::abs(m)) > l )|| l > l_max_)
  {
    return;
  }
  if (m >= 0)
  {
    values_[m][l] =  x;
  }
  else
  {
    values_[m][l] =  is_even(m) ? x : -x;
  }
}

void NormalizedLegendreArray::unsafe_set(const natural_t l, const natural_t m, const real_t x)
{
  values_[m][l] =  x;
}

NormalizedLegendreArray::NormalizedLegendreArray(NormalizedLegendreArray &&other) noexcept :
  l_max_(other.l_max_)
{
  values_ = std::move(other.values_);
}

NormalizedLegendreArray &NormalizedLegendreArray::operator=(NormalizedLegendreArray &&other) noexcept
{
  if (this != &other)
  {
    l_max_ = other.l_max_;
    values_ = std::move(other.values_);
  }
  return *this;
}

natural_t NormalizedLegendreArray::l_max() const 
{
  return l_max_;
}

}

