#include <iostream>
#include <vector>
#include <complex>
#include <fftw3.h>
#include <omp.h>
#include <chrono>
#include <atomic>
#include <algorithm>
#include <ipp.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <experimental/unordered_map>

//#include <stdio.h>

#define EXIT_MAIN exitLine:                                  
#define check_sts(st) if((st) != ippStsNoErr) goto exitLine;
 
using namespace std;

const double PI = 3.1415926535897932384;

int main() {
    // Sampling frequency
    const double fs = 100000;
    // Sampling interval
    const double T = 1.0 / fs;
    // Wave frequency
    const double f = 100;
   // const double f2 = 600;
    // Number of samples
    const int n = 1024*256;
    //fftw_init_threads();
    //fftw_plan_with_nthreads(1);
    // Allocate memory for the input and output arrays
    fftw_complex* in, * out;
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
    double* indeneme = (double*)malloc(sizeof(double) * n);
    

    
    tbb::cache_aligned_allocator<double> deneme;
        
    double* RealIn;
    double* RealOut;
    RealIn  = (double*)fftw_alloc_real(sizeof(double)*n);
    RealOut = (double*)fftw_alloc_real(sizeof(double)*n);
   // omp_allocator_handle_t On_dram = omp_low_lat_mem_alloc;
    
   // void* deneme = omp_alloc(sizeof(double) * n, On_dram);

    //lock->_lk = omp_alloc(sizeof(double) * n, On_dram);
    //omp_init_lock(lock);
    // double* HBW_MeM_Start= (double*)omp_alloc(sizeof(double) * n, omp_high_bw_mem_alloc);
    //fftw_r2r_kind RealKind;
    fftw_plan plan = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan RealPlan = fftw_plan_r2r_1d(n, RealIn, RealOut, FFTW_HC2R, FFTW_ESTIMATE);
    double frequency;
    double rms;
    int l = 0;
    int k = 0;
    int m = 0;
    double t;
    int index = 0;
    vector<double> magnitude(n);
    Ipp64f* magDst = ippsMalloc_64f(n * sizeof(Ipp64f));
    Ipp64f* SourceDst = ippsMalloc_64f(n * sizeof(Ipp64f));
    double* mag = (double*)malloc(sizeof(double) * n);
    // double magni[n];
    // Generate the input signal
    int i;

    IppStatus status;
       #pragma omp parallel for  
        for (i = 0; i < n; i++) {
             t = i * T;
            in[i][0] = sin(2 * PI * f * t);// +sin(2 * PI * (f2)+t);
            //indeneme[i] = sin(2 * PI * f * t);
           // RealIn[i] = sin(2 * PI * f * t);
             // RealIn++;
           in[i][1] = 0;
        }
        
        // Create the FFT plan
        
        //#pragma omp parallel for 
        auto start = chrono::high_resolution_clock::now();
       for (l = 0; l < 150; l++)
        {
            // Execute the FFT 
            fftw_execute(plan);
            
           #pragma omp parallel for
            
            for (int i = 0; i < n; i++)
            {
                SourceDst[i] = (out[i][0] * out[i][0] + out[i][1] * out[i][1]) ;
            }
            double deneme = SourceDst[0];
            check_sts(status = ippsSqrt_64f(SourceDst,magDst, n));
            EXIT_MAIN
            #pragma omp parallel for
                for (int i = 0; i < n; i++)
                {
                    magDst[i] = magDst[i] / n;
                }
            //for (k = 0; k < n; k++) {
            //  //  mag[k] = //sqrt(out[k][0] * out[k][0] + out[k][1] * out[k][1]) / n;
            //    
            //    
            //}
            
           //  Find the index of the maximum magnitude
            double max = magDst[0];
            #pragma omp parallel for //shared(index,magnitude) // private (m)
            for (m = 1; m < n; m++) {
                if (magDst[m] > max) {
                    max = magDst[m];
                    index = m;
                }
            }
             auto it = std::minmax_element(magnitude.begin(), magnitude.end());
             index = std::distance(magnitude.begin(), it.second);

            //  Calculate the frequency

            frequency = index / (n * T);

         //    Calculate the RMS value
            double sum = 0.0;
           #pragma omp parallel for  reduction(+:sum) 
            for (i = 0; i < n; i++) {
                //  #pragma omp atomic 
                sum += in[i][0] * in[i][0];
            }
            rms = sqrt(sum / n);   
    }
    // Print the results
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    std::cout << "Time taken by function: " << duration.count() << "millisec" << endl;
    std::cout << "Frequency: " << frequency << " Hz" << std::endl;
    std::cout << "RMS value: " << rms << std::endl;
    // Clean up
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
    return 0;
}
