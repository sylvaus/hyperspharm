#include <gtest/gtest.h>
#include <gsl/gsl_sf.h>
#include <algorithm>
#include "gegenbauer.h"

using hyperspharm::natural_t;
using hyperspharm::integer_t;
using hyperspharm::real_t;
using hyperspharm::Factorial;
using hyperspharm::GegenbauerPoly;
using hyperspharm::GegenbauerArray;

namespace
{

class GegenbauerTest : public ::testing::Test {
protected:
  GegenbauerTest()
  {
    for (real_t i = X_MIN; i <= X_MAX; i += (X_MAX - X_MIN)/(NB_X - 1))
    {
      x_values.push_back(i);
    }

    for (natural_t l = 0; l <= L_MAX; ++l)
    {
      l_values.push_back(l);
    }

    for (integer_t m = 0; m <= M_MAX; ++m)
    {
      m_values.push_back(m);
    }
  }


public:
  static constexpr real_t X_MIN = -1;
  static constexpr real_t X_MAX = 1;
  static constexpr natural_t NB_X = 100;

  static constexpr natural_t L_MAX = 500;
  static constexpr natural_t NB_L = 501;
  static constexpr integer_t M_MAX = 500;
  static constexpr natural_t NB_M = 501;

  std::vector<real_t> x_values;
  std::vector<natural_t> l_values;
  std::vector<integer_t> m_values;

  static real_t normalization_factor(const natural_t l, const natural_t m)
  {
    // sqrt(2 ** (2*m-1) * (l + m) / pi)
    auto first_part = std::sqrt((std::pow(2.0, (2.0 * m - 1.0)) * (m + l)) /  M_PI);

    // l! * ((m-1)! ** 2) / (2 * m + l - 1)
    std::vector<real_t> denominators;
    std::vector<real_t> nominators;
    for (natural_t i = 1; i <= (m - 1); ++i)
    {
      nominators.push_back(i);
      nominators.push_back(i);
    }

    for (natural_t i = l+1; i <= (2 * m) + l - 1; ++i)
    {
      denominators.push_back(i);
    }

    real_t second_part = 1;
    for (size_t j = 0; j < std::min(denominators.size(), nominators.size()); ++j)
    {
      second_part *= (static_cast<real_t>(nominators[j])/static_cast<real_t>(denominators[j]));
    }

    for (size_t j = std::min(denominators.size(), nominators.size()); j < nominators.size(); ++j)
    {
      second_part *= static_cast<real_t>(nominators[j]);
    }

    for (size_t j = std::min(denominators.size(), nominators.size()); j < denominators.size(); ++j)
    {
      second_part /= static_cast<real_t>(denominators[j]);
    }

    return first_part * sqrt(second_part);
  }
};

TEST_F(GegenbauerTest, SingleValueLowOrder)
{
  for (auto x : x_values)
  {
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(0, 1, x) * normalization_factor(0, 1), GegenbauerPoly::get_normalized(0, 1, x));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(1, 1, x) * normalization_factor(1, 1), GegenbauerPoly::get_normalized(1, 1, x));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(1, 2, x) * normalization_factor(1, 2), GegenbauerPoly::get_normalized(1, 2, x));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(0, 2, x) * normalization_factor(0, 2), GegenbauerPoly::get_normalized(0, 2, x));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(0, 10, x) * normalization_factor(0, 10), GegenbauerPoly::get_normalized(0, 10, x));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(1, 10, x) * normalization_factor(1, 10), GegenbauerPoly::get_normalized(1, 10, x));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(2, 2, x) * normalization_factor(2, 2), GegenbauerPoly::get_normalized(2, 2, x));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(5, 5, x) * normalization_factor(5, 5), GegenbauerPoly::get_normalized(5, 5, x));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(10, 2, x) * normalization_factor(10, 2), GegenbauerPoly::get_normalized(10, 2, x));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(10, 10, x) * normalization_factor(10, 10), GegenbauerPoly::get_normalized(10, 10, x));
  }
}

TEST_F(GegenbauerTest, SingleValueHighOrder)
{
  for (auto x : x_values)
  {
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(0, 50, x) * normalization_factor(0, 50), GegenbauerPoly::get_normalized(0, 50, x));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(25, 50, x) * normalization_factor(25, 50), GegenbauerPoly::get_normalized(25, 50, x));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(17, 85, x) * normalization_factor(17, 85), GegenbauerPoly::get_normalized(17, 85, x));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(0, 115, x) * normalization_factor(0, 115), GegenbauerPoly::get_normalized(0, 115, x));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(67, 115, x) * normalization_factor(67, 115), GegenbauerPoly::get_normalized(67, 115, x));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(112, 98, x) * normalization_factor(112, 98), GegenbauerPoly::get_normalized(112, 98, x));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(45, 87, x) * normalization_factor(45, 87), GegenbauerPoly::get_normalized(45, 87, x));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(98, 45, x) * normalization_factor(98, 45), GegenbauerPoly::get_normalized(98, 45, x));
  }
}

TEST_F(GegenbauerTest, GetArrayLowOrder)
{
  for (auto x : x_values)
  {
    auto array = GegenbauerPoly::get_norm_array(10, x, 1);
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(0, 1, x) * normalization_factor(0, 1), array.get(0, 1));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(1, 1, x) * normalization_factor(1, 1), array.get(1, 1));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(1, 2, x) * normalization_factor(1, 2), array.get(1, 2));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(0, 2, x) * normalization_factor(0, 2), array.get(0, 2));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(0, 10, x) * normalization_factor(0, 10), array.get(0, 10));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(1, 10, x) * normalization_factor(1, 10), array.get(1, 10));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(2, 2, x) * normalization_factor(2, 2), array.get(2, 2));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(5, 5, x) * normalization_factor(5, 5), array.get(5, 5));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(10, 2, x) * normalization_factor(10, 2), array.get(10, 2));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(10, 1, x) * normalization_factor(10, 1), array.get(10, 1));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(10, 10, x) * normalization_factor(10, 10), array.get(10, 10));
  }
}

TEST_F(GegenbauerTest, GetArrayHighOrder)
{
  for (auto x : x_values)
  {
    auto array = GegenbauerPoly::get_norm_array(115, x, 1);
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(0, 50, x) * normalization_factor(0, 50), array.get(0, 50));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(25, 50, x) * normalization_factor(25, 50), array.get(25, 50));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(17, 85, x) * normalization_factor(17, 85), array.get(17, 85));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(0, 115, x) * normalization_factor(0, 115), array.get(0, 115));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(67, 115, x) * normalization_factor(67, 115), array.get(67, 115));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(112, 98, x) * normalization_factor(112, 98), array.get(112, 98));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(45, 87, x) * normalization_factor(45, 87), array.get(45, 87));
    ASSERT_FLOAT_EQ(gsl_sf_gegenpoly_n(98, 45, x) * normalization_factor(98, 45), array.get(98, 45));
  }
}

}
