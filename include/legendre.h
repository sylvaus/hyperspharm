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

class NormalizedLegendreArray;

class LegendrePoly
{
public:
  static const real_t SQRT_4_PI;
  static const real_t SPHARM_NORM;
  static const real_t FULLY_NORM;

  
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

  static NormalizedLegendreArray get_norm_array(const real_t normalization_coeff,
                                                const natural_t l_max,
                                                const real_t x);

  static NormalizedLegendreArray get_fully_norm_array(const natural_t l_max,
                                                      const real_t x);
  static NormalizedLegendreArray get_sph_norm_array(const natural_t l_max,
                                                    const real_t x);
  
private:

  static void check_parameters(const natural_t l, 
                               const integer_t m, 
                               const real_t x);
};

class NormalizedLegendreArray
{
friend LegendrePoly;

public:
  const natural_t l_max;

  explicit NormalizedLegendreArray(const natural_t l_max);

  real_t get(const natural_t l, const integer_t m);
  void set(const natural_t l, const integer_t m, const real_t x);

private:
  size_t get_index(const natural_t l, const integer_t m);
  std::vector<real_t> values;
};

}
