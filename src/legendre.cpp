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

const real_t LegendrePoly::sqrt_4_pi = sqrt(4.0 * M_PI);

real_t LegendrePoly::get(const natural_t l, const real_t x)
{
  return get_associated(l, 0, x);
}


real_t LegendrePoly::get_associated(const natural_t l, 
                                    const integer_t m, 
                                    const real_t x)
{  
#ifndef NOCHECK
  check_parameters(l, m, x);
#endif
  if (almost_equal(std::abs(x), static_cast<real_t>(1.0), static_cast<real_t>(0.000001))) 
  {
    if (m != 0) {return 0;}
    if (x > 0.0) {return 1.0;}
    else {return minus_one_power(l);}
  }
 
  const real_t sqrt_1_x2 = sqrt(1.0d - x*x);
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
    p_l_1_m = ((k-1) * (inv_half_sqrt_1_x2) * p_l_1_m) - ((l - (k-2)) * (l + (k - 1)) * p_l__m);
    p_l__m = tmp;
  }

  if (m == abs_m)
  {
    p_l_1_m = minus_one_power(abs_m) * (Factorial::get(l + abs_m)/Factorial::get(l - abs_m)) * p_l_1_m;
  }

  return (minus_one_power(m) == 1.0) ? p_l_1_m : -p_l_1_m;
}


real_t LegendrePoly::get_fully_normalized(const natural_t l, 
                                          const integer_t m, 
                                          const real_t x)
{
#ifndef NOCHECK
  check_parameters(l, m, x);
#endif
  const auto abs_m = static_cast<natural_t >(std::abs(m));
  
  // Compute N_m^m
  real_t n_m_m = 1;
  if (m != 0)
  {
    const real_t poly = (1.0 - (x * x));
    for (natural_t i = 1; i <= abs_m; i++)
    {
      n_m_m *= (poly * 
                (2.0 * static_cast<real_t>(i) + 1.0) /
                 (2.0 * static_cast<real_t>(i)));
    }
    n_m_m = minus_one_power(abs_m) * std::sqrt(n_m_m);
  }
  
  if (l == abs_m) {return n_m_m;}
  if (l == -abs_m) 
  {
    if (is_even(abs_m)) {return n_m_m;}
    else {return -n_m_m;}
  }
  
  // Compute N_{m + 1}^m
  const auto abs_m_real = static_cast<real_t>(abs_m);
  
  real_t n_m_1_m = x * sqrt((2.0 * abs_m_real) + 3.0) * n_m_m;
  for (natural_t i = (abs_m + 2); i <= l; i++)
  {
    const real_t coeff_n_m_1_m = sqrt(static_cast<real_t>((4 * i * i) - 1) / // (2*i - 1) * (2*i + 1) = (2*i)²  - 1 = (4 * i²) -1
                                      static_cast<real_t>((i - abs_m) * (i + abs_m))
    );
    const real_t coeff_n_m_m = sqrt(static_cast<real_t>(((2 * i) + 1) * (i - abs_m - 1) * (i + abs_m - 1))/
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
  return get_fully_normalized(l, m, x) / sqrt_4_pi;
}


void LegendrePoly::check_parameters(const natural_t l, 
                                      const integer_t m, 
                                      const real_t x)
{
  const auto abs_m = static_cast<natural_t>(std::abs(m));
  if (abs_m > l)
  {
    throw std::invalid_argument( "Associated Legendre Polynomial: abs(m) must be smaller or equal to l" );
  }

  if (std::abs(x) > 1.0)
  {
    throw std::invalid_argument( "Associated Legendre Polynomial: function domain is x in [-1, 1] " );
  }
}

}

