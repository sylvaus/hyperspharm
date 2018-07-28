#include "spharms.h"
#include "gtest/gtest.h"

using namespace hyperspharm;

TEST(SphericalSurface, InitValueConstructor)
{
  const unsigned int n = 1024;
  const unsigned int m = 2048;
  SphericalSurface surface(n, m, 4.0);
  for (natural_t i = 0; i < surface.rows(); ++i)
  {
    for (natural_t j = 0; j < surface.cols(); ++j)
    {
      EXPECT_FLOAT_EQ(surface.get(i, j), 4.0);
    }
  }
}

TEST(SphericalSurface, MapFuncRealVoid)
{
  const unsigned int n = 1024;
  const unsigned int m = 2048;
  SphericalSurface surface(n, m);
  std::function<real_t()> func = []() { return 4.0; };
  surface.map(func);
  for (natural_t i = 0; i < surface.rows(); ++i)
  {
    for (natural_t j = 0; j < surface.cols(); ++j)
    {
      EXPECT_FLOAT_EQ(surface.get(i, j), 4.0);
    }
  }
}

TEST(SphericalSurface, MapFuncRealOldValue)
{
  const unsigned int n = 1024;
  const unsigned int m = 2048;
  SphericalSurface surface(n, m, 2.0);
  std::function<real_t(const real_t)> func =
      [](const real_t old_val) { return old_val * 2.0; };
  surface.map(func);

  for (natural_t i = 0; i < surface.rows(); ++i)
  {
    for (natural_t j = 0; j < surface.cols(); ++j)
    {
      EXPECT_FLOAT_EQ(surface.get(i, j), 4.0);
    }
  }
}

TEST(SphericalSurface, MapFuncRealNMOldValue)
{
  const unsigned int n = 1024;
  const unsigned int m = 2048;
  SphericalSurface surface(n, m, 2.0);
  std::function<real_t(const natural_t, const natural_t, const real_t)> func =
      [](const natural_t theta_n, const natural_t psi_m, const real_t old_val) {
          return old_val * (theta_n * psi_m);
      };
  surface.map(func);

  for (natural_t theta_n = 0; theta_n < surface.rows(); ++theta_n)
  {
    for (natural_t psi_m = 0; psi_m < surface.cols(); ++psi_m)
    {
      EXPECT_FLOAT_EQ(surface.get(theta_n, psi_m), 2.0 * (theta_n * psi_m));
    }
  }
}