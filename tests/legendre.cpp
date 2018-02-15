#include "legendre.h"
#include "gsl/gsl_sf.h"
#include "gtest/gtest.h"
#include <chrono>
#include <iostream>

namespace hyperspharm
{

/**
 * @brief Evaluate Legendre Polynomial 2,1 for the value x
 * This function is used to check if the Legendre Polynomial implementation
 * follows the Wikipedia notation: https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
 * 
 * @param x value for which the Legendre Polynomial 2,1 should be evaluated
 * @return hyperspharm::real_t
 */
real_t P_2_1(real_t x)
{
  return -3.0 * x * pow(1.0 - pow(x, 2), 0.5);
}
  
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
  EXPECT_FLOAT_EQ(LegendrePoly::get(3, 0), 0);
  EXPECT_FLOAT_EQ(LegendrePoly::get(3, 1), 1);
  EXPECT_FLOAT_EQ(LegendrePoly::get(4, 0), 3.0/8.0);
  EXPECT_FLOAT_EQ(LegendrePoly::get(4, 1), 1.0);
  EXPECT_FLOAT_EQ(LegendrePoly::get(10, 0), -0.24609375);
  EXPECT_FLOAT_EQ(LegendrePoly::get(10, 0.5), -0.1882286072);
  EXPECT_FLOAT_EQ(LegendrePoly::get(10, 1), 1);
  EXPECT_FLOAT_EQ(LegendrePoly::get(20, 1), 1);
  EXPECT_FLOAT_EQ(LegendrePoly::get(30, 1), 1);
  EXPECT_FLOAT_EQ(LegendrePoly::get(50, 1), 1);
  EXPECT_FLOAT_EQ(LegendrePoly::get(100, 1), 1);
}

TEST(AssociatedLegendre, NotationCheck) 
{
  EXPECT_FLOAT_EQ(P_2_1(-0.5), LegendrePoly::get_associated(2, 1, -0.5));
  EXPECT_FLOAT_EQ(P_2_1(0.5), LegendrePoly::get_associated(2, 1, 0.5));  
}

TEST(AssociatedLegendre, LowOrder) 
{
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(2, 1, -1), gsl_sf_legendre_Plm(2, 1, -1));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(2, 1, -0.5), gsl_sf_legendre_Plm(2, 1, -0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(2, 1, 0), gsl_sf_legendre_Plm(2, 1, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(2, 1, 0.5), gsl_sf_legendre_Plm(2, 1, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(2, 1, 1), gsl_sf_legendre_Plm(2, 1, 1));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(3, 2, -1), gsl_sf_legendre_Plm(3, 2, -1));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(3, 2, -0.5), gsl_sf_legendre_Plm(3, 2, -0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(3, 2, 0), gsl_sf_legendre_Plm(3, 2, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(3, 2, 0.5), gsl_sf_legendre_Plm(3, 2, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(3, 2, 1), gsl_sf_legendre_Plm(3, 2, 1));
}

TEST(AssociatedLegendre, HighOrder) 
{
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(100, 20, -0.5), gsl_sf_legendre_Plm(100, 20, -0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(100, 20, 0.5), gsl_sf_legendre_Plm(100, 20, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(100, 19, -0.5), gsl_sf_legendre_Plm(100, 19, -0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(100, 19, 0.5), gsl_sf_legendre_Plm(100, 19, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(120, 40, -0.9), gsl_sf_legendre_Plm(120, 40, -0.9));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(120, 40, 0.3), gsl_sf_legendre_Plm(120, 40, 0.3));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(120, 49, -0.7), gsl_sf_legendre_Plm(120, 49, -0.7));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(120, 49, 0.2), gsl_sf_legendre_Plm(120, 49, 0.2));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(90, 80, -0.8), gsl_sf_legendre_Plm(90, 80, -0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(90, 80, 0.5), gsl_sf_legendre_Plm(90, 80, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(90, 89, -0.5), gsl_sf_legendre_Plm(90, 89, -0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_associated(90, 89, 0.5), gsl_sf_legendre_Plm(90, 89, 0.5));
}

}
