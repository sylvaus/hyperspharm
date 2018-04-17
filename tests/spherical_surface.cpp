#include "spharms.h"
#include "gtest/gtest.h"

using namespace hyperspharm;

TEST(SphericalSurface, InitValueConstructor)
{
  const unsigned int size = 1024;
  SphericalSurface surface(size, size, 4.0);
  for (natural_t i = 0; i < size; ++i)
  {
    for (natural_t j = 0; j < size; ++j)
    {
      EXPECT_FLOAT_EQ(surface.get(i, j), 4.0);
    }
  }
}

TEST(SphericalSurface, MapFuncRealVoid)
{
  const unsigned int size = 1024;
  SphericalSurface surface(size, size);
  surface.map([](){return 4.0;});
  for (natural_t i = 0; i < size; ++i)
  {
    for (natural_t j = 0; j < size; ++j)
    {
      EXPECT_FLOAT_EQ(surface.get(i, j), 4.0);
    }
  }
}

TEST(SphericalSurface, MapFuncRealOldValue)
{
  const unsigned int size = 1024;
  SphericalSurface surface(size, size, 2.0);
  std::function<real_t(const real_t)> func =
      [](const real_t old_val){return old_val * 2.0;};
  surface.map(func);

  for (natural_t i = 0; i < size; ++i)
  {
    for (natural_t j = 0; j < size; ++j)
    {
      EXPECT_FLOAT_EQ(surface.get(i, j), 4.0);
    }
  }
}

TEST(SphericalSurface, MapFuncRealIJOldValue)
{
  const unsigned int size = 1024;
  SphericalSurface surface(size, size, 2.0);
  std::function<real_t(const natural_t, const natural_t, const real_t)> func =
      [](const natural_t theta_n, const natural_t psi_m, const real_t old_val)
      {
        return old_val * (theta_n * psi_m);
      };
  surface.map(func);

  for (natural_t theta_n = 0; theta_n < size; ++theta_n)
  {
    for (natural_t psi_m = 0; psi_m < size; ++psi_m)
    {
      EXPECT_FLOAT_EQ(surface.get(theta_n, psi_m), 2.0 * (theta_n * psi_m));
    }
  }
}