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
  EXPECT_FLOAT_EQ(Factorial::get(200),78865786736479050355236321393218506229513597768717326329474253324435944996340334292030428401198462390417721213891963883025764279024263710506192662495282993111346285727076331723739698894392244562145166.0);
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
  EXPECT_FLOAT_EQ(Factorial::get(200),78865786736479050355236321393218506229513597768717326329474253324435944996340334292030428401198462390417721213891963883025764279024263710506192662495282993111346285727076331723739698894392244562145166.0);
}

TEST(PrimeFactors, SmallNonPrime) 
{
  natural_t number = 12;
  std::vector<natural_t> expected_result = {2, 2, 3};
  auto result = PrimeFactors::compute(number);
  EXPECT_EQ(result.size(), expected_result.size());
  for (unsigned int index = 0; index < expected_result.size(); index++)
  {
    EXPECT_EQ(result[index], expected_result[index]);
  }
}

TEST(PrimeFactors, SmallPrime) 
{
  natural_t number = 13;
  std::vector<natural_t> expected_result = {13};
  auto result = PrimeFactors::compute(number);
  EXPECT_EQ(result.size(), expected_result.size());
  for (unsigned int index = 0; index < expected_result.size(); index++)
  {
    EXPECT_EQ(result[index], expected_result[index]);
  }
}

TEST(PrimeFactors, BigNonPrime) 
{
  natural_t number = 1245487;
  std::vector<natural_t> expected_result = {31, 40177};
  auto result = PrimeFactors::compute(number);
  EXPECT_EQ(result.size(), expected_result.size());
  for (unsigned int index = 0; index < expected_result.size(); index++)
  {
    EXPECT_EQ(result[index], expected_result[index]);
  }
}

TEST(PrimeFactors, BigPrime) 
{
  natural_t number = 104729;
  std::vector<natural_t> expected_result = {104729};
  auto result = PrimeFactors::compute(number);
  EXPECT_EQ(result.size(), expected_result.size());
  for (unsigned int index = 0; index < expected_result.size(); index++)
  {
    EXPECT_EQ(result[index], expected_result[index]);
  }
}

}
