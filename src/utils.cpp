/**
 * @file utils.cpp
 * @author Sylvaus
 * @date Sat Jan 13 2018
 * @brief Utilities functions and classes
 */

#include "utils.h"

namespace hyperspharm
{
  
const natural_t Factorial::MAX_MEMOISATION_N = 1024;
std::vector<real_t> Factorial::factorials_ = {0};
  
real_t Factorial::get(const natural_t n)
{
  if (factorials_.size() > n)
  {
  return factorials_[n];
  }
  
  for (natural_t index = factorials_.size(); index <= n; index++)
  {
    factorials_.push_back(factorials_[index - 1] * static_cast<real_t>(index));
  }
  
  return factorials_[n];
}

const natural_t PrimeFactors::MAX_MEMOISATION_N = 1024;
natural_t PrimeFactors::max_prime_checked_ = 12;
std::vector<natural_t> PrimeFactors::primes_ = {2, 3, 5, 7};

std::vector<natural_t> PrimeFactors::compute(const natural_t n)
{
  natural_t sqrt_n = std::floor(std::sqrt(n));
  
  bool max_memo_reached = false;
  if ((max_prime_checked_ < sqrt_n))
  {
     max_memo_reached = compute_primes(sqrt_n);
  }
  
  std::vector<natural_t> divisors;
  if (max_memo_reached)
  {
    // TODO: Implement factorization for number bigger than MAX_MEMOISATION_N ** 2
  }
  else
  {
    for (auto divisor : primes_)
    {
      if(sqrt_n < divisor)
      {
        break;
      }
      if ((n % divisor) == 0)
      {
        // TODO: Implement finding order of divisor and divisor bigger tha sqrt(n)
        divisors.push_back(divisor);
      }
    }
  }
  
  
  return divisors;
}

bool PrimeFactors::compute_primes(const natural_t n)
{
  if (primes_.size() >= MAX_MEMOISATION_N)
  {
    return true;
  }
  // TODO: Can be greatly optimized 
  bool max_memo_reached = false;
  for (natural_t i = max_prime_checked_ + 1; i <= n; i+= 2)
  {
    std::cout << i << std::endl;
    max_prime_checked_ = i;
    if (is_prime(i))
    {
      primes_.push_back(i);
      if (primes_.size() >= MAX_MEMOISATION_N)
      {
        max_memo_reached = true;
        break;
      }
    }
  }
  
  max_prime_checked_++;
  return max_memo_reached;
}

bool PrimeFactors::is_prime(const natural_t n)
{
  // This function assumes that all the primes number smaller or equal to
  // sqrt(n) are known
  natural_t sqrt_n = std::floor(std::sqrt(n));
  for (auto divisor : primes_)
  {
    if(sqrt_n < divisor)
    {
      return true;
    }
    if ((n % divisor) == 0)
    {
      return false;
    }
  }
  return true;
}

  
}
