/**
 * @file utils.h
 * @author Sylvaus
 * @date Sat Jan 13 2018
 * @brief Utilities functions and classes
 */

#pragma once

#include <vector>
#include <map>
#include <cmath>
#include <array>
#include <iostream>

#include "types.h"

namespace hyperspharm
{

template<class T>
bool almost_equal(T x, T y, T tolerance)
{
    return std::abs(x-y) <= tolerance;
}
 
class Factorial
{
public:
  static real_t get(const natural_t n);
  
private:
  static const natural_t MAX_MEMOISATION_N;
  static std::vector<real_t> factorials_;
};

class PrimeFactors
{
public: 
  static std::vector<natural_t> compute(natural_t n);
  static bool is_prime(const natural_t n);

private:
  static const natural_t MAX_MEMOISATION_N;
  static natural_t max_prime_checked_;
  
  static std::vector<natural_t> primes_;
  
  static bool is_prime_memo(const natural_t n);
  static bool compute_primes(const natural_t n);
  inline static void check_divisor(natural_t &n, 
                                   const natural_t divisor, 
                                   std::vector<natural_t> &divisors);
};

}
