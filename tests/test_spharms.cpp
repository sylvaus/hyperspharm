#include "spharms.h"
#include "gtest/gtest.h"

using namespace hyperspharm;

TEST(Spharms, Order0)
{
  const unsigned int size = 256;
  SphericalSurface surface(size, size, 1.0);
  auto result = Spharm::spharm_transform(surface);
  EXPECT_FLOAT_EQ(result.get(0, 0).real(), 2.0 * std::sqrt(M_PI));
  EXPECT_FLOAT_EQ(result.get(0, 0).imag(), 0.0);

  for (natural_t l = 1; l < result.l_max() ; ++l)
  {
    for (natural_t m = 0; m <= l; ++m)
    {

    }
  }
}

