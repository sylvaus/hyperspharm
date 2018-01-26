/**
 * @file legendre.cpp
 * @author Sylvaus
 * @date Mon Jan 08 2018
 * @brief 
 *
 * Legendre Polynomial Helpers
 */

#include "legendre.h"

namespace hyperspharm
{
  
std::vector<std::vector<real_t>> LegendrePoly::alk_ = {{1}, {1}};  


real_t LegendrePoly::get(const natural_t l, const real_t value)
{
  compute_all_alk(l);
  
  real_t k = ((l % 2) == 0) ? 0 : 1;
  real_t result = 0;
  for (real_t coeff : alk_[l])
  {
    result += coeff * pow(value, k);
    k += 2.0;
  }
  return result;
}


real_t LegendrePoly::get_associated(const natural_t l, 
                                    const integer_t m, 
                                    const real_t value)
{
  real_t result = 0;
  const natural_t abs_m = std::abs(m);
  
  for (natural_t k = 0; k <= (l-abs_m); k++)
  {
    result += (Factorial::get(k + abs_m) / Factorial::get(k)) * coeff(l, m+k) * pow(value, k);
  }
  
  if (m < 0)
  {
    if ((m % 2) == 0) 
    {
      result *= (Factorial::get(l + m) / Factorial::get(l - m));
    }
    else
    {
      result *= -(Factorial::get(l + m) / Factorial::get(l - m));
    }
  }
  
  return pow((1 - value * value), static_cast<real_t>(abs_m) / 2.0) * result;
}

real_t LegendrePoly::coeff(const natural_t l, const natural_t k)
{
  compute_all_alk(l);
  if ((l % 2) == 0)
  {
    return ((k % 2) == 0) ? alk_[l/2][0] : 0;
  }
  else
  {
    return ((k % 2) == 1) ? alk_[(l-1)/2][0] : 0;
  }
}

/**
 * @brief Compute the polynomial coefficients of all Legendre polynomial of order inferior to l_max
 * 
 * The Legendre Polynomial can be expressed with the following formula:
 *   P_l(x) = sum_{k = 0}^{l} a_{l,k} x^k
 * where  a_{l,k} =
 *   * 0 if l+k is odd
 *   * (-1)^((l-k)/2) * ((l+k)(l+k-1)...(k+1))/(2^l * l!) 
 *         * binomial coefficient (l, (l+k)/2)
 * 
 * Since we know the coefficients for n = 0, 1, we can use Bonnet’s recursion formula to find the next coefficients
 *   (l + 1)P_{l + 1}(x) = (2l + 1)xP_l(x) - lP_{l-1}(x)
 *   => P_{l}(x) = (2l-1) xP_{l-1}(x) / l - (l-1)P_{l-2}(x) / l
 *   => * a_{l, 0} = -(l-1)a_{l-2, 0} / l
 *      * a_{l, l} = (2l-1)a_{l-1, l-1} / l
 *      * a_{l, l-1} = 0 since l + l-1 is odd
 *      * a_{l, k} = (2l-1)a_{l-1, k-1} / l -(l-1)a_{l-2, k} / l for the other k
 * 
 * Moreover, only the non zero values are saved
 * 
 * TODO: Add limitation on the size of alk_
 * 
 * @param l_max maximum order of Legendre Polynomial 
 */
void LegendrePoly::compute_all_alk(const natural_t l_max)
{
  for (natural_t l = alk_.size(); l <= l_max; l++)
  {
    const real_t p_l_2_coeff = -(static_cast<real_t>(l) - 1.0) / static_cast<real_t>(l);
    const real_t p_l_1_coeff = (2.0 * static_cast<real_t>(l) - 1.0) / static_cast<real_t>(l);
    
    std::vector<real_t> ak;
    if ((l % 2) == 0)
    {
      const natural_t nb_non_zero_coeffs = (l/2) + 1;
      ak.reserve(nb_non_zero_coeffs);
      ak.push_back(p_l_2_coeff * alk_[l-2][0]);
      for (natural_t k = 1; k < nb_non_zero_coeffs - 1; k++)
      {
        ak.push_back((p_l_1_coeff * alk_[l-1][k-1]) + (p_l_2_coeff * alk_[l-2][k])); 
      }
      ak.push_back(p_l_1_coeff * alk_[l-1][nb_non_zero_coeffs - 2]);
    }
    else
    {
      const natural_t nb_non_zero_coeffs = (l+1) / 2;
      ak.reserve(nb_non_zero_coeffs);
      for (natural_t k = 0; k < nb_non_zero_coeffs - 1; k++)
      {
        ak.push_back((p_l_1_coeff * alk_[l-1][k]) + (p_l_2_coeff * alk_[l-2][k])); 
      }
      ak.push_back(p_l_1_coeff * alk_[l-1][nb_non_zero_coeffs - 1]);
    }
    
    alk_.push_back(ak);
  }
}

}

