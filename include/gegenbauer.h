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

  real_t get(natural_t l, integer_t m) const;
  real_t unsafe_get(natural_t l, natural_t m) const;
  void set(natural_t l,integer_t m, real_t x);
  void unsafe_set(natural_t l, natural_t m, real_t x);

  natural_t l_max() const;
private:
  natural_t l_max_;
  std::vector<std::vector<real_t>> values_;
};



class GegenbauerPoly
{
public:
  static real_t get(natural_t l, integer_t m, real_t x);
  static GegenbauerArray get_array(natural_t l, real_t x);

};

}



