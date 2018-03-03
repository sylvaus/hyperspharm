/**
 * @file Benchmark.cpp
 * @author Sylvaus
 * @date March 03 2018
 * @brief Compare implementation speed against other implementations
 */

#include <chrono>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <random>
#include <gsl/gsl_sf.h>
#include "legendre.h"


struct legendre_test_value
{
  hyperspharm::natural_t l;
  hyperspharm::natural_t m;
  hyperspharm::real_t x;
};

std::vector<legendre_test_value> generate_legendre_test_values(const unsigned long nb_values, const hyperspharm::natural_t max_l)
{
  std::vector<legendre_test_value> test_values;
  test_values.reserve(nb_values);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> x_dis(0.0, 1.0);
  std::uniform_int_distribution<> l_dis(0, max_l);

  for(unsigned long i = 0;  i < nb_values; i++)
  {
    legendre_test_value value;
    value.l = static_cast<hyperspharm::natural_t>(l_dis(gen));
    value.m = static_cast<hyperspharm::natural_t>(l_dis(gen));
    value.m = (value.l == 0) ? 0 : (value.m % value.l);
    value.x = x_dis(gen);
    test_values.push_back(value);
  }

  return test_values;
}

void test_spharm_normalized_legendre()
{
  const unsigned int nb_tests = 10000;
  const hyperspharm::natural_t max_l = 10000;
  std::vector<legendre_test_value> values = generate_legendre_test_values(nb_tests, max_l);
  std::vector<hyperspharm::real_t> results_this, results_gnu;
  results_this.reserve(nb_tests);
  results_gnu.reserve(nb_tests);

  auto start_this = std::chrono::high_resolution_clock::now();
  for(const auto& value : values)
  {
    results_this.push_back(hyperspharm::LegendrePoly::get_spharm_normalized(value.l,value.m, value.x));
  }
  auto finish_this = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_time_this = finish_this - start_this;

  auto start_gnu = std::chrono::high_resolution_clock::now();
  for(const auto& value : values)
  {
    results_gnu.push_back(gsl_sf_legendre_sphPlm(value.l, value.m, value.x));
  }
  auto finish_gnu = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_time_gnu = finish_gnu - start_gnu;

  std::cout << "gsl_sf_legendre_sphPlm took " << elapsed_time_gnu.count()
            << "s to complete " << nb_tests << " computations \n";

  std::cout << "LegendrePoly::get_spharm_normalized took " << elapsed_time_this.count()
            << "s to complete " << nb_tests << " computations \n";
}

int main()
{
  test_spharm_normalized_legendre();
  return 0;
}