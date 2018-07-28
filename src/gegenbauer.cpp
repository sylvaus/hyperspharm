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

namespace hyperspharm
{

real_t GegenbauerPoly::get(natural_t , integer_t , real_t )
{
  return 0;
}

}

