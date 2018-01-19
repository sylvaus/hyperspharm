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
std::vector<real_t> Factorial::factorials_ = {1};
  
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
natural_t PrimeFactors::max_prime_checked_ = 7;
std::vector<natural_t> PrimeFactors::primes_ = {2, 3, 5, 7};

std::vector<natural_t> PrimeFactors::compute(natural_t n)
{
  const natural_t sqrt_n = std::floor(std::sqrt(n));
  
  bool max_memo_reached = false;
  if ((max_prime_checked_ < sqrt_n))
  {
     max_memo_reached = compute_primes(sqrt_n);
  }
  
  std::vector<natural_t> divisors;
  for (auto divisor : primes_)
  {
    if ((sqrt_n < divisor) || (n < divisor)) {break;}
    check_divisor(n, divisor, divisors);
  }
  
  if (max_memo_reached)
  {
    std::cout << "Careful: PrimeFactors::compute was not designed for large number \n"
              << "         the computation may take a really long time\n";
    for (natural_t divisor = max_prime_checked_ + 1; divisor < sqrt_n; divisor++)
    {
      if ((sqrt_n < divisor) || (n < divisor)) {break;}
      if (is_prime(divisor))
      {
        check_divisor(n, divisor, divisors);
      }
    }
  }
  
  if (n > 1) {divisors.push_back(n);}
  
  return divisors;
}

bool PrimeFactors::is_prime(const natural_t n)
{
  if (n < 2) {return false;}
  if (2 == n) {return true;}
  natural_t sqrt_n = std::floor(std::sqrt(n));
  for (natural_t divisor = 3; divisor < sqrt_n; divisor += 2)
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

bool PrimeFactors::is_prime_memo(const natural_t n)
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

bool PrimeFactors::compute_primes(const natural_t n)
{
  static natural_t wheel_index = 0;
  const std::array<natural_t, 8> wheel_index_values = {{4, 2, 4, 2, 4, 6, 2, 6}};
  
  if (primes_.size() >= MAX_MEMOISATION_N)
  {
    return true;
  }
  while (max_prime_checked_ < n)
  {
    max_prime_checked_ += wheel_index_values[wheel_index];
    if ((++wheel_index) == 8) {wheel_index = 0;}
    
    std::cout << max_prime_checked_ << std::endl;
    if (is_prime_memo(max_prime_checked_))
    {
      primes_.push_back(max_prime_checked_);
      if (primes_.size() >= MAX_MEMOISATION_N)
      {
        return true;
        break;
      }
    }
  }
  
  return false;
}

inline void PrimeFactors::check_divisor(natural_t &n,
                                        const natural_t divisor, 
                                        std::vector<natural_t> &divisors)
{
  auto div = std::ldiv(n, divisor);
  while (div.rem == 0)
  {
    n = div.quot;
    divisors.push_back(divisor);
    div = std::ldiv(n, divisor);
  }
}

}
