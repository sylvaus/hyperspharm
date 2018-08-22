#include <gegenbauer.h>
#include <gtest/gtest.h>

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

  static real_t normalization_factor(const natural_t l, const integer_t m)
  {
    return std::sqrt((std::pow(2, (2 * m - 1)) * (m + l) * Factorial::get(l)) /  (M_PI * Factorial::get(2 * m + l - 1)) ) 
            * Factorial::get(m - 1);
  }

  static real_t Gegenbauer_l0_mAny()
  {
    return 1;
  }

  static real_t Gegenbauer_l1(integer_t m, real_t x)
  {
    return 2 * m * x;
  }
};

TEST_F(GegenbauerTest, SingleValueLowOrder)
{
  for (auto x : x_values)
  {
    ASSERT_FLOAT_EQ(Gegenbauer_l0_mAny() * normalization_factor(0, 1), GegenbauerPoly::get_normalized(0, 1, x));
    std::cout << x; 
  }
  
  for (auto x : x_values)
  {
    ASSERT_FLOAT_EQ(Gegenbauer_l1(1, x) * normalization_factor(1, 1), GegenbauerPoly::get_normalized(1, 1, x));
    std::cout << x;
  }

}

}
