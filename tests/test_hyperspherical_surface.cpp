#include "hyperspharm.h"
#include "gtest/gtest.h"

using namespace hyperspharm;

TEST(HyperSphericalSurface, InitValueConstructor)
{
  const natural_t n = 92;
  const natural_t l = 120;
  const natural_t m = 152;
  HyperSphericalSurface surface(n, l, m, 4.0);
  for (natural_t i = 0; i < surface.theta_nb(); ++i)
  {
    for (natural_t j = 0; j < surface.psi_nb(); ++j)
    {
      for (natural_t k = 0; k < surface.phi_nb(); ++k)
      {
        EXPECT_FLOAT_EQ(surface.get(i, j, k), 4.0);
      }
    }
  }
}

TEST(HyperSphericalSurface, MapFuncRealVoid)
{
  const natural_t n = 92;
  const natural_t l = 120;
  const natural_t m = 152;
  HyperSphericalSurface surface(n, l, m);
  std::function<real_t()> func = []() { return 4.0; };
  surface.map(func);
  for (natural_t i = 0; i < surface.theta_nb(); ++i)
  {
    for (natural_t j = 0; j < surface.psi_nb(); ++j)
    {
      for (natural_t k = 0; k < surface.phi_nb(); ++k)
      {
        EXPECT_FLOAT_EQ(surface.get(i, j, k), 4.0);
      }
    }
  }
}

TEST(HyperSphericalSurface, MapFuncRealOldValue)
{
  const natural_t n = 92;
  const natural_t l = 120;
  const natural_t m = 152;
  HyperSphericalSurface surface(n, l, m, 2.0);
  std::function<real_t(const real_t)> func =
      [](const real_t old_val) { return old_val * 2.0; };
  surface.map(func);

  for (natural_t i = 0; i < surface.theta_nb(); ++i)
  {
    for (natural_t j = 0; j < surface.psi_nb(); ++j)
    {
      for (natural_t k = 0; k < surface.phi_nb(); ++k)
      {
        EXPECT_FLOAT_EQ(surface.get(i, j, k), 4.0);
      }
    }
  }
}

TEST(HyperSphericalSurface, MapFuncRealNLMOldValue)
{
  const natural_t n = 92;
  const natural_t l = 120;
  const natural_t m = 152;
  HyperSphericalSurface surface(n, l, m, 2.0);
  std::function<real_t(const natural_t, const natural_t, const natural_t, const real_t)> func =
      [](const natural_t theta_i, const natural_t psi_i, const natural_t phi_i, const real_t old_val) {
        return old_val * (theta_i * std::sqrt(psi_i) / (phi_i + 1));
      };
  surface.map(func);

  for (natural_t i = 0; i < surface.theta_nb(); ++i)
  {
    for (natural_t j = 0; j < surface.psi_nb(); ++j)
    {
      for (natural_t k = 0; k < surface.phi_nb(); ++k)
      {
        EXPECT_FLOAT_EQ(surface.get(i, j, k), 2.0 * (i * std::sqrt(j) / (k + 1)));
      }
    }
  }
}

