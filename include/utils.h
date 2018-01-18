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

#include "types.h"

namespace hyperspharm
{
 
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
  static std::vector<natural_t> compute(const natural_t n);

private:
  static const natural_t MAX_MEMOISATION_N;
  static natural_t max_prime_checked_;
  
  static std::vector<natural_t> primes_;
  
  static bool compute_primes(const natural_t n);
  static bool is_prime(const natural_t n);
};

}
