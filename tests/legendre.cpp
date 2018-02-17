#include "legendre.h"
#include "gsl/gsl_sf.h"
#include "gtest/gtest.h"
#include <chrono>
#include <iostream>

using hyperspharm::LegendrePoly;

namespace hyperspharm_legendrepoly_test
{
const hyperspharm::real_t sqrt_4_pi = sqrt(4.0 * M_PI);
  
/**
 * @brief Evaluate Legendre Polynomial 2,1 for the value x
 * This function is used to check if the Legendre Polynomial implementation
 * follows the Wikipedia notation: https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
 * 
 * @param x value for which the Legendre Polynomial 2,1 should be evaluated
 * @return hyperspharm::real_t
 */
hyperspharm::real_t P_2_1(hyperspharm::real_t x)
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

TEST(FullyNormalizedAssociatedLegendre, LowOrderSectorial) // l==m 
{
  
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(1, 1, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(1, 1, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(1, 1, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(1, 1, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(1, 1, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(1, 1, 1));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(2, 2, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(2, 2, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(2, 2, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(2, 2, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(2, 2, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(2, 2, 1));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(3, 3, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(3, 3, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(3, 3, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(3, 3, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(3, 3, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(3, 3, 1));
}

TEST(FullyNormalizedAssociatedLegendre, LowOrderNonSectorial) // l!=m 
{
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(1, 0, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(1, 0, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(1, 0, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(1, 0, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(1, 0, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(1, 0, 1));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(2, 1, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(2, 1, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(2, 1, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(2, 1, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(2, 1, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(2, 1, 1));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(2, 0, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(2, 0, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(2, 0, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(2, 0, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(2, 0, 0.8) / sqrt_4_pi, gsl_sf_legendre_sphPlm(2, 0, 0.8));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(2, 0, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(2, 0, 1));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(3, 2, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(3, 2, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(3, 2, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(3, 2, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(3, 2, 0.8) / sqrt_4_pi, gsl_sf_legendre_sphPlm(3, 2, 0.8));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(3, 2, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(3, 2, 1));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(3, 1, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(3, 1, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(3, 1, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(3, 1, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(3, 1, 0.8) / sqrt_4_pi, gsl_sf_legendre_sphPlm(3, 1, 0.8));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(3, 1, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(3, 1, 1));
}

TEST(FullyNormalizedAssociatedLegendre, HighOrderSectorial) // l==m 
{
  
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(100, 100, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(100, 100, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(100, 100, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(100, 100, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(100, 100, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(100, 100, 1));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(200, 200, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(200, 200, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(200, 200, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(200, 200, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(200, 200, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(200, 200, 1));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(300, 300, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(300, 300, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(300, 300, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(300, 300, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(300, 300, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(300, 300, 1));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(101, 101, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(101, 101, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(101, 101, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(101, 101, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(101, 101, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(101, 101, 1));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(201, 201, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(201, 201, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(201, 201, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(201, 201, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(201, 201, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(201, 201, 1));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(301, 301, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(301, 301, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(301, 301, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(301, 301, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(301, 301, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(301, 301, 1));
}

TEST(FullyNormalizedAssociatedLegendre, HighOrderNonSectorial) // l!=m 
{
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(100, 0, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(100, 0, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(100, 0, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(100, 0, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(100, 0, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(100, 0, 1));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(200, 100, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(200, 100, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(200, 100, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(200, 100, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(200, 100, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(200, 100, 1));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(200, 0, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(200, 0, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(200, 0, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(200, 0, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(200, 0, 0.8) / sqrt_4_pi, gsl_sf_legendre_sphPlm(200, 0, 0.8));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(200, 0, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(200, 0, 1));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(300, 200, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(300, 200, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(300, 200, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(300, 200, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(300, 200, 0.8) / sqrt_4_pi, gsl_sf_legendre_sphPlm(300, 200, 0.8));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(300, 200, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(300, 200, 1));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(300, 100, 0) / sqrt_4_pi, gsl_sf_legendre_sphPlm(300, 100, 0));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(300, 100, 0.5) / sqrt_4_pi, gsl_sf_legendre_sphPlm(300, 100, 0.5));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(300, 100, 0.8) / sqrt_4_pi, gsl_sf_legendre_sphPlm(300, 100, 0.8));
  EXPECT_FLOAT_EQ(LegendrePoly::get_fully_normalized(300, 100, 1) / sqrt_4_pi, gsl_sf_legendre_sphPlm(300, 100, 1));
}

}
