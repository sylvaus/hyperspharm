
#include <chrono>
#include <iostream>
#include "fft.h"
#include "utils.h"
#include "types.h"

void test () {
    const hyperspharm::natural_t nSamples = 2048;
    //real_t nSeconds = 1.0;                      // total time for sampling
    //real_t sampleRate = nSamples / nSeconds;    // n Hz = n / second 
    //real_t freqResolution = sampleRate / nSamples; // freq step in FFT result
    hyperspharm::complex_t x[nSamples];                // storage for sample data
    hyperspharm::complex_t X[nSamples];                // storage for FFT answer
    const int nFreqs = 5;
    hyperspharm::real_t freq[nFreqs] = { 2, 5, 11, 17, 29 }; // known freqs for testing
    
    // generate samples for testing
    for(hyperspharm::natural_t i=0; i<nSamples; i++) {
        x[i] = hyperspharm::complex_t(0.,0.);
        // sum several known sinusoids into x[]
        for(int j=0; j<nFreqs; j++)
            x[i] += sin( 2*M_PI*freq[j]*i/nSamples );
        X[i] = x[i];        // copy into X[] for FFT work & result
    }
    // compute fft for this data
    hyperspharm::fft(X,nSamples);
    
    /*printf("  n\tx[]\tX[]\tf\n");       // header line
    // loop to print values
    for(int i=0; i<nSamples; i++) {
        printf("% 3d\t%+.3f\t%+.3f\t%g\n",
            i, x[i].real(), abs(X[i]), i*freqResolution );
    }*/
}

int main() 
{
    //const long cycle = 1000;
    auto start = std::chrono::high_resolution_clock::now();
    /*for (long index = 0; index < cycle; index++)
        test ();*/
    auto results = hyperspharm::PrimeFactors::compute(1545548575124);
    for (auto result : results)
      std::cout << result  << "\n";
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<hyperspharm::real_t> elapsed = finish - start;
    std::cout << "Code ran for " << elapsed.count() << "ms\n";
    start = std::chrono::high_resolution_clock::now();
    /*for (long index = 0; index < cycle; index++)
        test ();*/
    results = hyperspharm::PrimeFactors::compute(1545548575124);
    for (auto result : results)
      std::cout << result  << "\n";
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    std::cout << "Code ran for " << elapsed.count() << "ms\n";
}
