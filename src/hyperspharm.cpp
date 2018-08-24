/**
 * @file hyperspharm.cpp
 * @author Sylvaus
 * @date Tue July 03 2018
 * @brief
 *
 */

#include "hyperspharm.h"

namespace hyperspharm
{

HyperSphericalSurface::HyperSphericalSurface(const natural_t theta_nb, const natural_t psi_nb, const natural_t phi_nb):
  theta_nb_(theta_nb), psi_nb_(psi_nb), phi_nb_(phi_nb), values_(theta_nb * psi_nb * phi_nb)
{

}


HyperSphericalSurface::HyperSphericalSurface(const natural_t theta_nb, const natural_t psi_nb,
                                             const natural_t phi_nb, const real_t init_val):
    theta_nb_(theta_nb), psi_nb_(psi_nb),
    phi_nb_(phi_nb), values_(theta_nb * psi_nb * phi_nb, init_val)
{

}

real_t HyperSphericalSurface::get(natural_t theta_i, natural_t psi_i, natural_t phi_i) const
{
  return values_[get_index(theta_i, psi_i, phi_i)];
}

void HyperSphericalSurface::set(natural_t theta_i, natural_t psi_i, natural_t phi_i, real_t radius)
{
  // TODO Check could be added (and then an unsafe method)
  values_[get_index(theta_i, psi_i, phi_i)] = radius;
}

natural_t HyperSphericalSurface::theta_nb() const
{
  return theta_nb_;
}

natural_t HyperSphericalSurface::psi_nb() const
{
  return psi_nb_;
}

natural_t HyperSphericalSurface::phi_nb() const
{
  return phi_nb_;
}

std::vector<real_t> HyperSphericalSurface::thetas() const
{
  return std::vector<real_t>();
}

std::vector<real_t> HyperSphericalSurface::psis() const
{
  return std::vector<real_t>();
}

std::vector<real_t> HyperSphericalSurface::phis() const
{
  return std::vector<real_t>();
}

void HyperSphericalSurface::map(std::function<real_t()> func)
{
  for(auto& value : values_)
  {
    value = func();
  }
}

void HyperSphericalSurface::map(std::function<real_t(real_t)> func)
{
  for(auto& value : values_)
  {
    value = func(value);
  }
}

void HyperSphericalSurface::map(std::function<real_t(natural_t, natural_t,
                                                     natural_t, real_t)> func)
{
  natural_t n = 0;
  natural_t l = 0;
  natural_t m = 0;
  for (natural_t i = 0; i < (theta_nb_ * psi_nb_ * phi_nb_); ++i)
  {
    if (phi_nb_ == m)
    {
      ++l;
      m = 0;
    }
    if (psi_nb_ == l)
    {
      ++n;
      l = 0;
    }
    values_[i] = func(n, l, m, values_[i]);
    ++m;
  }
}

std::string HyperSphericalSurface::to_string() const
{
  return "TODO";
}

std::vector<real_t> &HyperSphericalSurface::values()
{
  return values_;
}

HyperSphericalCoeffs::HyperSphericalCoeffs(natural_t n_max, natural_t l_max, natural_t m_max):
  n_max_(n_max), l_max_(l_max), m_max_(m_max), values_(n_max * l_max * m_max)
{
}

HyperSphericalCoeffs::HyperSphericalCoeffs(natural_t n_max, natural_t l_max,
                                           natural_t m_max, complex_t init_val):
    n_max_(n_max), l_max_(l_max), m_max_(m_max), values_(n_max * l_max * m_max, init_val)
{
}

complex_t HyperSphericalCoeffs::get(natural_t n, natural_t l, natural_t m) const
{
  return values_[get_index(n, l, m)];
}

void HyperSphericalCoeffs::set(natural_t n, natural_t l, natural_t m, complex_t value)
{
  values_[get_index(n, l, m)] = value;
}

natural_t HyperSphericalCoeffs::n_max() const
{
  return n_max_;
}

natural_t HyperSphericalCoeffs::l_max() const
{
  return l_max_;
}

natural_t HyperSphericalCoeffs::m_max() const
{
  return m_max_;
}

void HyperSphericalCoeffs::map(std::function<complex_t()> func)
{
  for(auto& value : values_)
  {
    value = func();
  }
}

void HyperSphericalCoeffs::map(std::function<complex_t(complex_t)> func)
{
  for(auto& value : values_)
  {
    value = func(value);
  }
}

void HyperSphericalCoeffs::map(std::function<complex_t(natural_t, natural_t,
                                                       natural_t, complex_t)> func)
{
  natural_t n = 0;
  natural_t l = 0;
  natural_t m = 0;
  for (natural_t i = 0; i < values_.size(); ++i)
  {
    if (m_max_ == m)
    {
      ++l;
      m = 0;
    }
    if (l_max_ == l)
    {
      ++n;
      l = 0;
    }
    values_[i] = func(n, l, m, values_[i]);
    ++m;
  }
}

std::string HyperSphericalCoeffs::to_string()
{
  return "TODO";
}

std::vector<complex_t>& HyperSphericalCoeffs::values()
{
  return values_;
}

}