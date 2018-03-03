/**
 * @file legendre.cpp
 * @author Sylvaus
 * @date Mon Jan 08 2018
 * @brief 
 *
 * Legendre Polynomial Helpers
 */

#pragma once

#include <cstdint>
#include <vector>
#include "utils.h"
#include "types.h"

namespace hyperspharm
{

class LegendrePoly
{
public:
  static const real_t sqrt_4_pi;
  
  static real_t get(const natural_t l, const real_t x);
  static real_t get_associated(const natural_t l, 
                               const integer_t m, 
                               const real_t x);
  static real_t get_fully_normalized(const natural_t l, 
                                     const integer_t m, 
                                     const real_t x);
  static real_t get_spharm_normalized(const natural_t l, 
                                      const integer_t m, 
                                      const real_t x);
  
private:
  static void check_parameters(const natural_t l, 
                               const integer_t m, 
                               const real_t x);
};

}
