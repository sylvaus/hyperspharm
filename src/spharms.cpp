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


SphericalSurface::SphericalSurface(const natural_t rows, const natural_t cols) :
  rows_(rows), cols_(cols), values_(rows * cols)
{
}

real_t SphericalSurface::get(natural_t theta_n, const natural_t psi_m) const
{
  return values_[cols_ * theta_n + psi_m];
}

void SphericalSurface::set(const natural_t theta_n, const natural_t psi_m, const real_t radius_nm)
{
  values_[cols_ * theta_n + psi_m] = radius_nm;
}

std::vector<complex_t> SphericalSurface::get_psi_array(const natural_t theta_n) const
{
  auto start_index = values_.begin() + (theta_n * cols_);
  return std::vector<complex_t>(start_index, start_index + cols_);
}

natural_t SphericalSurface::rows() const
{
  return rows_;
}

natural_t SphericalSurface::cols() const
{
  return cols_;
}

std::string SphericalSurface::to_string()
{
  std::stringstream sstream;
  sstream << std::scientific << std::setprecision(3) << std::left;
  sstream << "Theta\\Psi";
  for (unsigned int psi = 0; psi < cols_ ; ++psi)
  {
    sstream << std::setw(5) << std::right << (psi * 2 + 1) << "/" <<
               std::setw(3) << std::left << (cols_ * 2) << "pi";
  }
  sstream << std::endl;

  for (unsigned int theta = 0; theta < rows_ ; ++theta)
  {
    sstream << std::setw(3) << std::right << (theta * 2 + 1) << "/" <<
               std::setw(3) << std::left << (rows_ * 2) << "pi : ";
    for (unsigned int psi = 0; psi < cols_; ++psi)
    {
      sstream << std::setw(10) << get(theta, psi) << " ";
    }
    sstream << std::endl;
  }
  return sstream.str();
}


SphericalHarmonics::SphericalHarmonics(const natural_t l_max) :
  l_max_(l_max), values_(((l_max + 2)*(l_max + 1))/2)
{
}

complex_t SphericalHarmonics::get(const natural_t l, const natural_t m) const
{
  const size_t index = ((l+1) * l)/2 + m;
  if(index >= values_.size())
  {
    return {0, 0};
  }

  return values_[index];
}

void SphericalHarmonics::set(const natural_t l, const natural_t m, const complex_t value)
{
  const size_t index = ((l+1) * l)/2 + m;
  if(index >= values_.size())
  {
    return;
  }
  values_[index] = value;
}

natural_t SphericalHarmonics::l_max() const
{
  return l_max_;
}

std::string SphericalHarmonics::to_string()
{
  std::stringstream sstream;
  sstream << std::scientific << std::setprecision(2);
  sstream << " l\\m  ";
  for (unsigned int m = 0; m < l_max_ ; ++m)
  {
    sstream << std::setw(11) << m << std::setw(9) << " ";
  }
  sstream << std::endl;

  for (unsigned int l = 0; l < l_max_ ; ++l)
  {
    sstream << std::setw(3) << l << " : ";
    for (unsigned int m = 0; m <= l; ++m)
    {
      sstream << std::setw(18) << get(l, m) << " ";
    }
    sstream << std::endl;
  }
  return sstream.str();
}

SphericalHarmonics Spharm::spharm_transform(const SphericalSurface &spherical_surface)
{
  std::vector<std::vector<complex_t>> fm_thetas;
  fm_thetas.reserve(spherical_surface.rows());
  for (natural_t theta_index = 0; theta_index < spherical_surface.rows(); ++theta_index)
  {
    std::vector<complex_t> fm_theta = spherical_surface.get_psi_array(theta_index);
    fft(fm_theta.data(), fm_theta.size());
    fm_thetas.push_back(std::move(fm_theta));
  }

  std::vector<NormalizedLegendreArray> plm_thetas;
  plm_thetas.reserve(spherical_surface.rows());
  const real_t w = M_PI / static_cast<real_t>(spherical_surface.rows());
  const real_t delta_theta = M_PI / static_cast<real_t>(spherical_surface.rows());
  real_t theta = delta_theta / 2.0;
  for (natural_t theta_index = 0; theta_index < spherical_surface.rows(); ++theta_index)
  {
    NormalizedLegendreArray pnm_theta = LegendrePoly::get_sph_norm_array(spherical_surface.rows(), std::cos(theta));
    plm_thetas.push_back(std::move(pnm_theta));
    theta += delta_theta;
  }

  const auto fft_normalization = 2.0 * M_PI / static_cast<real_t>(spherical_surface.rows());
  SphericalHarmonics result(spherical_surface.rows());
  for (natural_t l = 0; l < result.l_max(); ++l)
  {
    for (natural_t m = 0; m <= l; ++m)
    {
      complex_t fm_n = {0, 0};
      theta = delta_theta / 2.0;
      for (natural_t theta_index = 0; theta_index < spherical_surface.rows(); ++theta_index)
      {
        fm_n += fm_thetas[theta_index][m] *
                (plm_thetas[theta_index].get(l, m) * w * std::sin(theta) * fft_normalization);
        theta += delta_theta;
      }
      result.set(l, m, fm_n);
    }
  }

  return result;
}
}
