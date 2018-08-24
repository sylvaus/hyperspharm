/**
 * @file hyperspharm.h
 * @author Sylvaus
 * @date Tue July 03 2018
 * @brief
 *
 */

#pragma once

#include <algorithm>
#include <vector>
#include <iomanip>
#include <cmath>
#include <functional>
#include "fft.h"
#include "utils.h"
#include "types.h"
#include "legendre.h"

namespace hyperspharm
{

class HyperSphericalSurface
{
public:
  HyperSphericalSurface(natural_t theta_nb, natural_t psi_nb, natural_t phi_nb);
  HyperSphericalSurface(natural_t theta_nb, natural_t psi_nb, natural_t phi_nb, real_t init_val);

  real_t get(natural_t theta_i, natural_t psi_i, natural_t phi_i) const;
  void set(natural_t theta_i, natural_t psi_i, natural_t phi_i, real_t radius);

  natural_t theta_nb() const;
  natural_t psi_nb() const;
  natural_t phi_nb() const;

  std::vector<real_t> thetas() const;
  std::vector<real_t> psis() const;
  std::vector<real_t> phis() const;

  void map(std::function<real_t ()>);
  void map(std::function<real_t (real_t old_val)>);
  void map(std::function<real_t (natural_t theta_i, natural_t psi_i, natural_t phi_i, real_t old_val)>);

  std::string to_string() const;

  inline natural_t get_index(natural_t theta_i, natural_t psi_i, natural_t phi_i) const
  {
    return (theta_i * (psi_nb_ * phi_nb_)) + (psi_i * phi_nb_) + phi_i;
  }
  std::vector<real_t>& values();
private:
  natural_t theta_nb_;
  natural_t psi_nb_;
  natural_t phi_nb_;
  std::vector<real_t> values_;
};

class HyperSphericalCoeffs
{
public:
  explicit HyperSphericalCoeffs(natural_t n_max, natural_t l_max, natural_t m_max);
  explicit HyperSphericalCoeffs(natural_t n_max, natural_t l_max,
                                natural_t m_max, complex_t init_val);

  // TODO: Add possibility to input m negative
  complex_t get(natural_t n, natural_t l, natural_t m) const;
  void set(natural_t n, natural_t l, natural_t m, complex_t value);

  natural_t n_max() const;
  natural_t l_max() const;
  natural_t m_max() const;

  void map(std::function<complex_t ()>);
  void map(std::function<complex_t (complex_t old_val)>);
  void map(std::function<complex_t (natural_t n, natural_t l, natural_t m, complex_t old_val)>);


  std::string to_string();

  inline natural_t get_index(natural_t n, natural_t l, natural_t m) const
  {
    return (n * (m_max_ * l_max_)) + (l * m_max_) + m;
  }
  std::vector<complex_t>& values();
private:
  natural_t n_max_;
  natural_t l_max_;
  natural_t m_max_;
  std::vector<complex_t> values_;
};

class HyperSpharm
{
  HyperSphericalCoeffs transform(HyperSphericalSurface surface);
  HyperSphericalSurface transform(HyperSphericalCoeffs coeffs);
};

}

