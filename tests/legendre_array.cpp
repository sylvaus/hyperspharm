#include "legendre.h"
#include "gsl/gsl_sf.h"
#include "gtest/gtest.h"
#include <chrono>
#include <iostream>

using hyperspharm::LegendrePoly;
using hyperspharm::NormalizedLegendreArray;

namespace hyperspharm_legendrepolyarray_test
{
const hyperspharm::real_t sqrt_4_pi = sqrt(4.0 * M_PI);


TEST(FullyNormalizedAssociatedLegendreArray, PositiveX) // l==m
{
  const hyperspharm::natural_t l_max = 1000;
  hyperspharm::real_t x = 0.5;
  auto *gsl_results = new double[gsl_sf_legendre_array_n(l_max)];
  gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_FULL, l_max,  x, -1, gsl_results);
  auto this_results = LegendrePoly::get_fully_norm_array(l_max, x);

  for (hyperspharm::natural_t l = 0; l <= l_max; ++l)
  {
    for (hyperspharm::natural_t m = 0; m <= l; ++m)
    {
      EXPECT_FLOAT_EQ(this_results.get(l, m), gsl_results[gsl_sf_legendre_array_index(l, m)]);
    }
  }

  delete[] gsl_results;
}

TEST(FullyNormalizedAssociatedLegendreArray, NegativeX) // l==m
{
  const hyperspharm::natural_t l_max = 1000;
  hyperspharm::real_t x = -0.5;
  auto *gsl_results = new double[gsl_sf_legendre_array_n(l_max)];
  gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_FULL, l_max,  x, -1, gsl_results);
  auto this_results = LegendrePoly::get_fully_norm_array(l_max, x);

  for (hyperspharm::natural_t l = 0; l <= l_max; ++l)
  {
    for (hyperspharm::natural_t m = 0; m <= l; ++m)
    {
      EXPECT_FLOAT_EQ(this_results.get(l, m), gsl_results[gsl_sf_legendre_array_index(l, m)]);
    }
  }

  delete[] gsl_results;
}


TEST(SpharmNormalizedAssociatedLegendreArray, PositiveX) // l!=m
{
  const hyperspharm::natural_t l_max = 1000;
  hyperspharm::real_t x = 0.5;
  auto *gsl_results = new double[gsl_sf_legendre_array_n(l_max)];
  gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM, l_max,  x, -1, gsl_results);
  auto this_results = LegendrePoly::get_sph_norm_array(l_max, x);

  for (hyperspharm::natural_t l = 0; l <= l_max; ++l)
  {
    for (hyperspharm::natural_t m = 0; m <= l; ++m)
    {
      EXPECT_FLOAT_EQ(this_results.get(l, m), gsl_results[gsl_sf_legendre_array_index(l, m)]);
    }
  }

  delete[] gsl_results;
}

TEST(SpharmNormalizedAssociatedLegendreArray, NegativeX) // l!=m
{
  const hyperspharm::natural_t l_max = 1000;
  hyperspharm::real_t x = 0.5;
  auto *gsl_results = new double[gsl_sf_legendre_array_n(l_max)];
  gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM, l_max,  x, -1, gsl_results);
  auto this_results = LegendrePoly::get_sph_norm_array(l_max, x);

  for (hyperspharm::natural_t l = 0; l <= l_max; ++l)
  {
    for (hyperspharm::natural_t m = 0; m <= l; ++m)
    {
      EXPECT_FLOAT_EQ(this_results.get(l, m), gsl_results[gsl_sf_legendre_array_index(l, m)]);
    }
  }

  delete[] gsl_results;
}
  
}
