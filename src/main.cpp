
#include <chrono>
#include <iostream>
#include "fft.h"

void test () {
    const uint32_t nSamples = 2048;
    //double nSeconds = 1.0;                      // total time for sampling
    //double sampleRate = nSamples / nSeconds;    // n Hz = n / second 
    //double freqResolution = sampleRate / nSamples; // freq step in FFT result
    std::complex<double> x[nSamples];                // storage for sample data
    std::complex<double> X[nSamples];                // storage for FFT answer
    const int nFreqs = 5;
    double freq[nFreqs] = { 2, 5, 11, 17, 29 }; // known freqs for testing
    
    // generate samples for testing
    for(uint32_t i=0; i<nSamples; i++) {
        x[i] = std::complex<double>(0.,0.);
        // sum several known sinusoids into x[]
        for(int j=0; j<nFreqs; j++)
            x[i] += sin( 2*M_PI*freq[j]*i/nSamples );
        X[i] = x[i];        // copy into X[] for FFT work & result
    }
    // compute fft for this data
    fft2(X,nSamples);
    
    /*printf("  n\tx[]\tX[]\tf\n");       // header line
    // loop to print values
    for(int i=0; i<nSamples; i++) {
        printf("% 3d\t%+.3f\t%+.3f\t%g\n",
            i, x[i].real(), abs(X[i]), i*freqResolution );
    }*/
}

int main() 
{
    const long cycle = 1000;
    auto start = std::chrono::high_resolution_clock::now();
    for (long index = 0; index < cycle; index++)
        test ();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Code ran for " << elapsed.count() << "ms\n";
}
