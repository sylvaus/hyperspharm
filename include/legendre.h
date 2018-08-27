/**
 * @file legendre.cpp
 * @author Sylvaus
 * @date Mon Jan 08 2018
 * @brief 
 *
 * Legendre Polynomial Helpers
 */

#pragma once

#include <cmath>
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
  static void compute_coefficients(const natural_t l_max);

  typedef struct {real_t alm; real_t blm;} coeff;
  static std::vector<std::vector<LegendrePoly::coeff>> coeffs_;
};

class NormalizedLegendreArray
{
friend LegendrePoly;

public:
  explicit NormalizedLegendreArray(const natural_t l_max);

  NormalizedLegendreArray (NormalizedLegendreArray&& other) noexcept;
  NormalizedLegendreArray& operator= (NormalizedLegendreArray&& other) noexcept;

  real_t get(const natural_t l, const integer_t m) const;
  real_t unsafe_get(const natural_t l, const natural_t m) const;
  void set(const natural_t l, const integer_t m, const real_t x);
  void unsafe_set(const natural_t l, const natural_t m, const real_t x);

  natural_t l_max() const;
private:
  natural_t l_max_;
  std::vector<std::vector<real_t>> values_;
};

}
