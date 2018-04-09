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

};

class NormalizedLegendreArray
{
friend LegendrePoly;

public:
  explicit NormalizedLegendreArray(const natural_t l_max);

  NormalizedLegendreArray (NormalizedLegendreArray&& other) noexcept;
  NormalizedLegendreArray& operator= (NormalizedLegendreArray&& other) noexcept;

  real_t get(const natural_t l, const integer_t m) const;
  void set(const natural_t l, const integer_t m, const real_t x);

  natural_t l_max();
private:
  static size_t get_index(const natural_t l, const integer_t m);
  natural_t l_max_;
  std::vector<real_t> values_;
};

}
