/**
 * @file types.h
 * @author Sylvaus
 * @date Mon Jan 08 2018
 * @brief File containing the definition of the different types
 *
 */

#pragma once 

#include <complex>
#include <vector>

namespace hyperspharm
{
  typedef uint64_t natural_t;
  typedef int64_t integer_t;
  typedef double real_t;
  typedef std::complex<real_t> complex_t;
  
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
#endif // !M_PI

}
