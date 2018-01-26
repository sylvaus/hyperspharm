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
  static real_t get(const natural_t l, const real_t value);
  static real_t get_associated(const natural_t l, 
                               const integer_t m, 
                               const real_t value);
  static real_t coeff(const natural_t l, const natural_t k);

private:
  static std::vector<std::vector<real_t>> alk_; // Legendre Polynomial Coefficents (only non zero-values)
  
  static void compute_all_alk(const natural_t l_max);
};

}
