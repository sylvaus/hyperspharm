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
    legendre_test_value value = {
        static_cast<hyperspharm::natural_t>(l_dis(gen))
        , static_cast<hyperspharm::natural_t>(l_dis(gen))
        , x_dis(gen)
    };
    value.m = (value.l == 0) ? 0 : (value.m % value.l);
    test_values.push_back(value);
  }

  return test_values;
}

void test_spharm_normalized_legendre()
{
  const unsigned int nb_tests = 1000000;
  const hyperspharm::natural_t max_l = 900;
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

  if (elapsed_time_gnu.count() > elapsed_time_this.count())
  {
    std::cout << ((elapsed_time_gnu.count() - elapsed_time_this.count()) / elapsed_time_gnu.count()) * 100
              << "% of improvement compared to reference\n";
  }
  else
  {
    std::cout << ((elapsed_time_this.count() - elapsed_time_gnu.count()) / elapsed_time_this.count()) * 100
              << "% worse compared to reference\n";
  }
}

void test_spharm_normalized_legendre_array()
{
  const unsigned int nb_tests = 1000;
  const hyperspharm::natural_t max_l = 900;
  std::vector<legendre_test_value> values = generate_legendre_test_values(nb_tests, max_l);
  std::vector<hyperspharm::real_t> results_this, results_gnu;
  results_this.reserve(nb_tests);
  results_gnu.reserve(nb_tests);

  auto start_this = std::chrono::high_resolution_clock::now();
  for(const auto& value : values)
  {
    auto this_results = hyperspharm::LegendrePoly::get_fully_norm_array(value.l, value.x);
    results_this.push_back(this_results.get(value.l,value.m));
  }
  auto finish_this = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_time_this = finish_this - start_this;

  auto start_gnu = std::chrono::high_resolution_clock::now();
  for(const auto& value : values)
  {
    auto *gsl_results = new double[gsl_sf_legendre_array_n(value.l)];
    gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_FULL, value.l,  value.x, -1, gsl_results);

    results_gnu.push_back(gsl_results[gsl_sf_legendre_array_index(value.l,value.m)]);
    delete[] gsl_results;
  }
  auto finish_gnu = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_time_gnu = finish_gnu - start_gnu;

  std::cout << "gsl_sf_legendre_sphPlm took " << elapsed_time_gnu.count()
            << "s to complete " << nb_tests << " computations \n";

  std::cout << "LegendrePoly::get_spharm_normalized took " << elapsed_time_this.count()
            << "s to complete " << nb_tests << " computations \n";
  if (elapsed_time_gnu.count() > elapsed_time_this.count())
  {
    std::cout << ((elapsed_time_gnu.count() - elapsed_time_this.count()) / elapsed_time_gnu.count()) * 100
              << "% of improvement compared to reference\n";
  }
  else
  {
    std::cout << ((elapsed_time_this.count() - elapsed_time_gnu.count()) / elapsed_time_this.count()) * 100
              << "% worse compared to reference\n";
  }
}

int main()
{
  test_spharm_normalized_legendre();
  test_spharm_normalized_legendre_array();
  return 0;
}