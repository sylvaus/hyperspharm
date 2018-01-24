/**
 * @file spharms.cpp
 * @author Sylvaus
 * @date Mon Jan 08 2018
 * @brief 
 *
 * Spherical Harmonics Helpers
 */

#include "spharms.h"

namespace hyperspharm
{
  
std::vector<std::vector<real_t>> LegendrePoly::alk_ = {{0}, {-0.5, 1.5}};  

/*
real_t LegendrePoly::get(const natural_t l, const complex_t value)
{
  return 0;
}

real_t LegendrePoly::get_associated(const natural_t l, 
                                    const natural_t m, 
                                    const complex_t value)
{
  return 0;
}*/

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
 *      * a_{l, k} = (2l-1)a_{l-1, k-1} / l - -(l-1)a_{l-2, k} / l for the other k
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
    ak.reserve(l+1);
    ak.push_back(p_l_2_coeff * alk_[l-2][0]);
    
    for (natural_t k = 1; k < l - 1; k++)
    {
      ak.push_back((p_l_1_coeff * alk_[l-1][k-1]) + (p_l_2_coeff * alk_[l-2][k])); 
    }
    ak.push_back(0);
    ak.push_back(p_l_1_coeff * alk_[l-1][l-1]);
    alk_.push_back(ak);
  }
  
}

}
