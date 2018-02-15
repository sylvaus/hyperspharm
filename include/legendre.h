/**
 * @file legendre.h
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

class LegendrePoly
{
public:
  static real_t get(const natural_t l, const real_t x);
  static real_t get_associated(const natural_t l, 
                               const integer_t m, 
                               const real_t x);
  static real_t get_fully_normalized(const natural_t l, 
                                     const integer_t m, 
                                     const real_t x);
  
private:
  static void check_parameters(const natural_t l, 
                               const integer_t m, 
                               const real_t x);
};

}
