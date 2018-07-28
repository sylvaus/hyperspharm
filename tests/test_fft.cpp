
#include "fft.h"
#include "gtest/gtest.h"
namespace hyperspharm
{
  
bool are_complex_equal(const std::complex<real_t>& left, const std::complex<real_t>& right)
{
  static const real_t threshold = 0.0001;
  return ((std::abs(left.real()- right.real()) < threshold) && (std::abs(left.imag()- right.imag()) < threshold));
}

TEST(FFT, FFTSizeNotPowerOfTwo) 
{
  const natural_t nSamples = 2048 + 1;
  std::complex<real_t> x[nSamples];
  EXPECT_FALSE(fft(x, nSamples));
}

TEST(FFT, IFFTSizeNotPowerOfTwo) 
{
  const natural_t nSamples = 2048 + 1;
  std::complex<real_t> x[nSamples];
  EXPECT_FALSE(ifft(x, nSamples));
}

TEST(FFT, ValidFFT) 
{
  const natural_t nSamples = 16;
  std::complex<real_t> x[nSamples] = {{1.0, 0}, {1.5, 0}, {-2.0, 0}, {2.5, 0}, {3.0, 0}, {-3.5, 0}, {4.0, 0}, {4.5, 0}, 
                                      {-5.0, 0}, {5.5, 0}, {6.0, 0}, {-6.5, 0}, {7.0, 0}, {7.5, 0}, {-8.0, 0}, {8.5, 0},};
  std::complex<real_t> expected_results[nSamples] = {{26, 0}, {-0.488467, 6.080799}, {0.142136, 1.899495}, {1.664545, - 0.202758}, {6, -2}, {38.619726, 13.454097}, {-28.142136, 17.899495}, {-15.795804, 3.737654}, {-14, 0}, {-15.795804, -3.737654}, {-28.142136, -17.899495}, {38.619726, -13.454097}, {6, 2}, {1.664545, 0.202758}, {0.142136, -1.899495}, {-0.488467, -6.080799}};
  EXPECT_TRUE(fft(x, nSamples));
  for(natural_t index = 0; index < nSamples; index++)
  {
    EXPECT_TRUE(are_complex_equal(x[index], expected_results[index]));
  }
}

TEST(FFT, ValidIFFT) 
{
  const natural_t nSamples = 16;
  std::complex<real_t> expected_results[nSamples] = {{1.0, 0}, {1.5, 0}, {-2.0, 0}, {2.5, 0}, {3.0, 0}, {-3.5, 0}, {4.0, 0}, {4.5, 0}, 
                                      {-5.0, 0}, {5.5, 0}, {6.0, 0}, {-6.5, 0}, {7.0, 0}, {7.5, 0}, {-8.0, 0}, {8.5, 0},};
  std::complex<real_t> x[nSamples] = {{26, 0}, {-0.488467, 6.080799}, {0.142136, 1.899495}, {1.664545, - 0.202758}, {6, -2}, {38.619726, 13.454097}, {-28.142136, 17.899495}, {-15.795804, 3.737654}, {-14, 0}, {-15.795804, -3.737654}, {-28.142136, -17.899495}, {38.619726, -13.454097}, {6, 2}, {1.664545, 0.202758}, {0.142136, -1.899495}, {-0.488467, -6.080799}};
  EXPECT_TRUE(ifft(x, nSamples));
  for(natural_t index = 0; index < nSamples; index++)
  {
    EXPECT_TRUE(are_complex_equal(x[index], expected_results[index]));
  }
}

}
