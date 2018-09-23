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

  /**
   * Get the value for order (l, m)
   * @param l
   * @param m
   * @return real_t
   * @throw invalid_argument if l or m is bigger than l_max
   */
  real_t get(natural_t l, natural_t m) const;
  /**
   * Same as get except that no check is done on (l, m).
   * The call may trigger a segfault if (l, m) is outside the domain [0..l_max]^2
   * @param l
   * @param m
   * @return
   */
  real_t unsafe_get(natural_t l, natural_t m) const;
  /**
   * Set the value for order (l, m)
   * @param l
   * @param m
   * @param x
   * @throw invalid_argument if l or m is bigger than l_max
   */
  void set(natural_t l,natural_t m, real_t x);
  /**
   * Same as set except that no check is done on (l, m).
   * The call may trigger a segfault if (l, m) is outside the domain [0..l_max]^2
   * @param l
   * @param m
   * @param x
   */
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
	static const real_t INV_SQRT_PI;
	static const real_t N01;

  /**
   * Returns the value of normalized G^m_l(x)
   * @param l natural number >= 0
   * @param m natural number > 0
   * @param x real number with domain -1 < x < 1
   * @return real
   * @throws invalid_argument if x, l, m outside their respective domains
   */
  static real_t get_normalized(natural_t l, natural_t m, real_t x);

  /**
   * Returns the values of normalized G^m_l(x) for l, m in [1, l_max]x[0, l_max]
   * @param l_max natural number: maximum order for the normalized Gegenbauer polynomial
   * @param x real number -1<x<1
   * @param coeff real number: all the values will multiplied by this coefficient
   * @return GegenbauerArray
   */
  static GegenbauerArray get_norm_array(natural_t l_max, real_t x, real_t coeff=1);
private:
  /**
   * Computes the recurrence coefficients and store them internally
   * @param l_max natural number: maximum order of the coefficients
   */
  static void compute_coefficients(natural_t l_max);

  typedef struct {real_t a; real_t b;} coeff;
  static std::vector<std::vector<GegenbauerPoly::coeff>> coeffs_;
};

}



