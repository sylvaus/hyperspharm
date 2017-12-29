/**
 * @file fft.cpp
 * @author Sylvaus
 * @date Thu Dec 28 2017
 * @brief fft library based on https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm#C++_Example_Code
 *
 */

#pragma once

#include <complex>
#include <stdio.h>

bool is_power_of_two(uint32_t value);
void separate (std::complex<double>* array, uint32_t size);
bool fft2 (std::complex<double>* array, uint32_t size);
void unsafe_fft2 (std::complex<double>* array, uint32_t size);
