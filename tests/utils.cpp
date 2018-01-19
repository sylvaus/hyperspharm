#include "utils.h"
#include "gtest/gtest.h"

namespace hyperspharm
{

TEST(Factorial, getWithoutMemo) 
{
  EXPECT_EQ(Factorial::get(0), 1);
  EXPECT_EQ(Factorial::get(1), 1);
  EXPECT_EQ(Factorial::get(2), 2);
  EXPECT_EQ(Factorial::get(3), 6);
  
  EXPECT_FLOAT_EQ(Factorial::get(22), 1124000727777607680000.0);
  EXPECT_FLOAT_EQ(Factorial::get(23), 25852016738884976640000.0);
}

TEST(Factorial, getWithMemo) 
{
  Factorial::get(23);
  Factorial::get(0);
  Factorial::get(1);
  Factorial::get(2);
  Factorial::get(3);
  
  EXPECT_EQ(Factorial::get(0), 1);
  EXPECT_EQ(Factorial::get(1), 1);
  EXPECT_EQ(Factorial::get(2), 2);
  EXPECT_EQ(Factorial::get(3), 6);
  
  EXPECT_FLOAT_EQ(Factorial::get(22), 1124000727777607680000.0);
  EXPECT_FLOAT_EQ(Factorial::get(23), 25852016738884976640000.0);
}
  
  
}
