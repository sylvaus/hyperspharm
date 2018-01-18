/**
 * @file fft.cpp
 * @author Sylvaus
 * @date Thu Dec 28 2017
 * @brief fft library based on https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm#C++_Example_Code
 *
 */

#pragma once

#include <iostream>
#include "types.h"

namespace hyperspharm
{

bool is_power_of_two(natural_t value);
void separate (complex_t* array, natural_t size);
bool fft (complex_t* array, natural_t size);
bool ifft (complex_t* array, natural_t size);
void unsafe_fft (complex_t* array, natural_t size);

}
