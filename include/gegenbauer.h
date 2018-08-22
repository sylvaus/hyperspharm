/**
 * @file gegenbauer.h
 * @author Sylvaus
 * @date Tue July 08 2018
 * @brief
 *      Gegenbauer Polynomial helper functions
 */

#pragma once

#include <cstdint>
#include <vector>
#include "utils.h"
#include "types.h"

namespace hyperspharm
{

class GegenbauerPoly;

class GegenbauerArray
{
  friend GegenbauerPoly;

public:
  explicit GegenbauerArray(const natural_t l_max);

  GegenbauerArray (GegenbauerArray&& other) noexcept;
  GegenbauerArray& operator= (GegenbauerArray&& other) noexcept;

  real_t get(natural_t l, natural_t m) const;
  real_t unsafe_get(natural_t l, natural_t m) const;
  void set(natural_t l,integer_t m, real_t x);
  void unsafe_set(natural_t l, natural_t m, real_t x);
  
  std::vector<std::vector<real_t>>& values();

  natural_t l_max() const;
private:
  natural_t l_max_;
  std::vector<std::vector<real_t>> values_;
};



class GegenbauerPoly
{
public:
  static constexpr real_t INV_SQRT_PI = 1.0 / std::sqrt(M_PI);
  static constexpr real_t N01 = std::sqrt(2.0) * INV_SQRT_PI;
  
  static real_t get_normalized(natural_t l, natural_t m, real_t x);
  static GegenbauerArray get_norm_array(natural_t l, real_t x);
private:
  static void compute_coefficients(const natural_t l_max);

  typedef struct {real_t alm; real_t blm;} coeff;
  static std::vector<std::vector<GegenbauerPoly::coeff>> coeffs_;
};

}



