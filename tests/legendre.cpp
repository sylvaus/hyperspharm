#include "legendre.h"
#include "gtest/gtest.h"

namespace hyperspharm
{

TEST(Legendre, InitialValues) 
{
  EXPECT_EQ(LegendrePoly::get(0, 0), 1);
  EXPECT_EQ(LegendrePoly::get(1, 0), 0);
  EXPECT_EQ(LegendrePoly::get(1, 1), 1);
}

TEST(Legendre, HigherOrder) 
{
  EXPECT_FLOAT_EQ(LegendrePoly::get(2, 0), -0.5);
  EXPECT_FLOAT_EQ(LegendrePoly::get(2, 1), 1);
  EXPECT_FLOAT_EQ(LegendrePoly::get(2, 2), 11.0/2.0);
  EXPECT_FLOAT_EQ(LegendrePoly::get(3, 0), 0);
  EXPECT_FLOAT_EQ(LegendrePoly::get(3, 1), 1);
  EXPECT_FLOAT_EQ(LegendrePoly::get(3, 2), 17);
  EXPECT_FLOAT_EQ(LegendrePoly::get(3, 3), 63);
  EXPECT_FLOAT_EQ(LegendrePoly::get(4, 0), 3.0/8.0);
  EXPECT_FLOAT_EQ(LegendrePoly::get(4, 1), 1.0);
  EXPECT_FLOAT_EQ(LegendrePoly::get(4, 3), 321);
  EXPECT_FLOAT_EQ(LegendrePoly::get(10, 0), -0.24609375);
  EXPECT_FLOAT_EQ(LegendrePoly::get(10, 0.5), -0.1882286072);
  EXPECT_FLOAT_EQ(LegendrePoly::get(10, 1), 1);
}

}
