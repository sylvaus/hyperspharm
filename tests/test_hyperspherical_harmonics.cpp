#include "hyperspharm.h"
#include "gtest/gtest.h"

using namespace hyperspharm;

TEST(HyperSphericalCoeffs, InitValueConstructor)
{
  const natural_t n = 92;
  const natural_t l = 120;
  const natural_t m = 152;
  HyperSphericalCoeffs coeffs(n, l, m, {4.0, 1.5});
  for (natural_t ni = 0; ni < coeffs.n_max(); ++ni)
  {
    for (natural_t li = 0; li < coeffs.l_max(); ++li)
    {
      for (natural_t mi = 0; mi < coeffs.m_max(); ++mi)
      {
        EXPECT_EQ(complex_t(4.0, 1.5), coeffs.get(ni, li, mi));
      }
    }
  }
}


TEST(HyperSphericalCoeffs, MapFuncRealVoid)
{
  const natural_t n = 3;
  const natural_t l = 2;
  const natural_t m = 2;
  HyperSphericalCoeffs coeffs(n, l, m);
  std::function<complex_t()> func = []() { return complex_t(3.14, 1.7); };
  coeffs.map(func);
  for (natural_t ni = 0; ni < coeffs.n_max(); ++ni)
  {
    for (natural_t li = 0; li < coeffs.l_max(); ++li)
    {
      for (natural_t mi = 0; mi < coeffs.m_max(); ++mi)
      {
        EXPECT_EQ(complex_t(3.14, 1.7), coeffs.get(ni, li, mi));
      }
    }
  }
}

TEST(HyperSphericalCoeffs, MapFuncRealOldValue)
{
  const natural_t n = 92;
  const natural_t l = 120;
  const natural_t m = 152;
  HyperSphericalCoeffs coeffs(n, l, m, complex_t(3.14, 1.7));
  std::function<complex_t(complex_t)> func =
      [](complex_t old_value) { return old_value * 2.0; };
  coeffs.map(func);
  for (natural_t ni = 0; ni < coeffs.n_max(); ++ni)
  {
    for (natural_t li = 0; li < coeffs.l_max(); ++li)
    {
      for (natural_t mi = 0; mi < coeffs.m_max(); ++mi)
      {
        EXPECT_EQ(complex_t(3.14 * 2.0, 1.7 * 2.0), coeffs.get(ni, li, mi));
      }
    }
  }
}

TEST(HyperSphericalCoeffs, MapFuncRealNLMOldValue)
{
  const natural_t n = 92;
  const natural_t l = 120;
  const natural_t m = 152;
  HyperSphericalCoeffs coeffs(n, l, m, complex_t(3.14, 1.7));
  std::function<complex_t(natural_t, natural_t, natural_t, complex_t)> func =
      [](natural_t ni, natural_t li, natural_t mi, complex_t old_val) {
        return old_val * (ni * std::sqrt(li) / (mi + 1));
      };
  coeffs.map(func);
  for (natural_t ni = 0; ni < coeffs.n_max(); ++ni)
  {
    for (natural_t li = 0; li < coeffs.l_max(); ++li)
    {
      for (natural_t mi = 0; mi < coeffs.m_max(); ++mi)
      {
        EXPECT_EQ(complex_t(3.14, 1.7) * (ni * std::sqrt(li) / (mi + 1)), coeffs.get(ni, li, mi));
      }
    }
  }
}
