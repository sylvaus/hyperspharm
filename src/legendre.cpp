/**
 * @file legendre.cpp
 * @author Sylvaus
 * @date Mon Jan 08 2018
 * @brief 
 *
 * Legendre Polynomial Helpers
 * Based on paper: http://www.scielo.org.co/pdf/racefn/v37n145/v37n145a09.pdf
 */

#include "legendre.h"

namespace hyperspharm
{
  
real_t minus_one_power(const natural_t m)
{
  return ((m % 2) == 0) ? 1.0 : -1.0;
}


real_t LegendrePoly::get(const natural_t l, const real_t x)
{
  return get_associated(l, 0, x);
}


real_t LegendrePoly::get_associated(const natural_t l, 
                                    const integer_t m, 
                                    const real_t x)
{
  if (almost_equal(std::abs(x), static_cast<real_t>(1.0), static_cast<real_t>(0.000001))) 
  {
    if (m != 0) {return 0;}
    if (x > 0.0) {return 1.0;}
    else {return minus_one_power(l);}
  }
  
  const natural_t abs_m = std::abs(m);
  if (abs_m > l)
  {
    throw std::invalid_argument( "Associated Legendre Polynomial: abs(m) must be smaller or equal to l" );
  }
  
  const real_t sqrt_1_x2 = sqrt(1.0 - x*x);
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
        return -0.5 * sqrt_1_x2;
      case 0:
        return x;
      case 1:
        return sqrt_1_x2;
    }
  }
  
  real_t p_l__m = minus_one_power(l) * pow(half_sqrt_1_x2, l) / Factorial::get(l); // P_l^{-l}
  if (m == -static_cast<integer_t>(l)) {return p_l__m;}
  real_t p_l_1_m = -p_l__m * static_cast<real_t>(l) * inv_half_sqrt_1_x2; // P_l^{1-l}
  
  // TODO: implement optimisation for m>0 using the relation between P_l^{m} and P_l^{-m} (https://en.wikipedia.org/wiki/Associated_Legendre_polynomials)
  for (integer_t k = 2 - l; k <= m; k++)
  {
    real_t tmp = p_l_1_m;
    p_l_1_m = ((k-1) * (inv_half_sqrt_1_x2) * p_l_1_m) - ((l - (k-2)) * (l + (k - 1)) * p_l__m);
    p_l__m = tmp;
  }
  return (minus_one_power(m) == 1.0) ? p_l_1_m : -p_l_1_m;
}

}

