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

}

