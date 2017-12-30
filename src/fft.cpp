/**
* @file fft.cpp
* @author Sylvaus
* @date Thu Dec 28 2017
* @brief fft library based on https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm#C++_Example_Code
*
*/

#include "fft.h"

using std::complex;

bool is_power_of_two(uint32_t value)
{
  while (((value % 2) == 0) && value > 1) 
  {
    value /= 2; 
  }
  return (value == 1);
}

/**
* @brief Separate even/odd elements to lower/upper halves of array respectively.
* 
* @param array array of complex number to be separated
* @param size size of the array
*/
void separate(std::complex< double >* array, uint32_t size)
{
  const uint32_t half_size = size / 2;
  complex<double>* temp_array = new complex<double>[half_size];  // get temp heap storage
  for(uint32_t index=0; index < half_size; index++)    // copy all odd elements to heap storage
    temp_array[index] = array[index * 2 + 1];
  for(uint32_t index=0; index < half_size; index++)    // copy all even elements to lower-half of a[]
    array[index] = array[index*2];
  for(uint32_t index=0; index < half_size; index++)    // copy all odd (from heap) to upper-half of a[]
    array[index+size/2] = temp_array[index];
  delete[] temp_array;                 // delete heap storage
}

/**
* @brief Compute the fft of the given array, only if size is a power of two.
* 
* Check if size is a power of two and then call unsafe_fft2.
* 
* @param array array of complex number to be separated
* @param size size of the array
* @return bool returns true if size is a power of two.
*/
bool fft (complex<double> array[], uint32_t size)
{
  if (is_power_of_two(size))
  {
    unsafe_fft(array, size);
    return true;
  }
  else
  {
    std::cerr << "Error in fft: Size of the array needs to be a power of 2\n"; 
    return false;
  }
}

/**
* @brief Compute the inverse fft of the given array, only if size is a power of two.
* 
* Check if size is a power of two and then call unsafe_fft2.
* 
* @param array array of complex number to be separated
* @param size size of the array
* @return bool returns true if size is a power of two.
*/
bool ifft (complex<double> array[], uint32_t size)
{
  if (is_power_of_two(size))
  {
    // Using the formula inverse F({X_n}) = F({X_{N-n}})/n;
    complex<double> temp;
    // X_n = X_{N-n}
    for (uint32_t index = 1; index < (size/2); index++)
    {
      temp = array[index];
      array[index] = array[size - index];
      array[size - index] = temp;
    }
    unsafe_fft(array, size);
    for (uint32_t index = 0; index < size; index++)
    {
      array[index] = array[index] / static_cast<double>(size);
    }
    return true;
  }
  else
  {
    std::cerr << "Error in ifft: Size of the array needs to be a power of 2\n"; 
    return false;
  }
}

/**
* @brief Compute the fft of the given array.
* 
* N must be a power-of-2, or bad things will happen.
* 
* N input samples in X[] are FFT'd and results left in X[].
* Because of Nyquist theorem, N samples means 
* only first N/2 FFT results in X[] are the answer.
* (upper half of X[] is a reflection with no new information).
* 
* @param array array of complex number to be separated
* @param size size of the array
*/
void unsafe_fft (complex<double> array[], uint32_t size) 
{
  if(size < 2) 
  {
      // bottom of recursion.
      // Do nothing here, because already X[0] = x[0]
  } 
  else 
  {
    const uint32_t half_size = size / 2;
    separate(array, size);      // all evens to lower half, all odds to upper half
    unsafe_fft(array, half_size);   // recurse even items
    unsafe_fft(array + half_size, half_size);   // recurse odd  items
    // combine results of two half recursions
    for(uint32_t index = 0; index < half_size; index++) 
    {
       complex<double> e = array[index    ];                              // even
       complex<double> o = array[index + half_size];                      // odd
       complex<double> w = exp( complex<double>(0,-2.*M_PI*index/size) ); // w is the "twiddle-factor"
       array[index] = e + w * o;
       array[index + half_size] = e - w * o;
    }
  }
}
