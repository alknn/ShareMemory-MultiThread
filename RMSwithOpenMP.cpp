//#include <iostream>
//#include <vector>
//#include <complex>
//#include <omp.h>
//#include <fftw3.h>
//#include <algorithm>
//#include <chrono>
//
//const double PI = 3.1415926535897932384;
//
//int main() {
//    // Sampling frequency
//    const double fs = 200000000;
//    // Sampling interval
//    const double T = 1.0 / fs;
//    // Wave frequency
//    const double f = 200000;
//    const double f2 = 210;
//    const double f4 = 510;
//    // Number of samples
//    const int n = 1024*128;
//    int num_peaks = 10;
//    double max_values[10];
//    int max_indices[10];
//    int N_THREADS = 8;
//    std::vector<double> magnitude(n);
//    fftw_init_threads();
//    fftw_plan_with_nthreads(8);
//
//    // Allocate memory for the input and output arrays
//    float* inputReal = fftwf_alloc_real(n);
//    fftw_complex* in, * out;
//    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
//    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
//    fftw_plan plan;
//    double frequency;
//    double rms;
//
//    // Generate the input signal
//    int i;
//    
//auto start = std::chrono::high_resolution_clock::now();
//#pragma omp parallel for shared(in,n)
//        for (i = 0; i < n; i++) {
//            double t = i * T;
//            in[i][0] = sin(2 * PI * f * t);
//            in[i][1] = 0;
//        }
//        
//        // Create the FFT plan
//         plan = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
//         fftw_make_planner_thread_safe();
//        
//        // Execute the FFT
//        fftw_execute(plan);
//     
//
//        // Calculate the magnitude
//       
//       
//#pragma omp parallel for shared(n)
//        for (i = 0; i < n; i++) {
//            magnitude[i] = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]) / n;
//        }
//        
//    
//        /*for (int i = 0; i < num_peaks; i++) {
//            max_values[i] = 0.0;
//            max_indices[i] = 0;
//        }
//        for (int i = 0; i < n / 2; i++) {
//            for (int j = 0; j < num_peaks; j++) {
//                if (magnitude[i] > max_values[j]) {
//                    for (int k = num_peaks - 1; k > j; k--) {
//                        max_values[k] = max_values[k - 1];
//                        max_indices[k] = max_indices[k - 1];
//                    }
//                    max_values[j] = magnitude[i];
//                    max_indices[j] = i;
//                    break;
//                }
//            }
//        }
//        for (int i = 0; i < num_peaks; i++) {
//            std::cout << (max_indices[i] * 100000 / n) << " Hz" << std::endl;
//        }*/
//        // Find the index of the maximum magnitude
//       
//        int index = 0;
//        double max = magnitude[0];
//#pragma omp parallel for shared (n)
//        for (i = 1; i < n; i++) {
//            if (magnitude[i] > max) {
//                max = magnitude[i];
//                index = i;
//            }
//        }
//        //std::sort(magnitude[0], magnitude[1027*64-1]);
//        // Calculate the frequency
//         frequency = index / (n * T);
//
//
//        // Calculate the RMS value
//        double sum = 0.0;
//#pragma omp parallel for reduction(+:sum)
//        for (i = 0; i < n; i++) {
//            sum += in[i][0] * in[i][0];
//        }
//         rms = sqrt(sum / n);
// 
//    auto end = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
//    std::cout << "Time taken by function: " << duration.count() << " millisec" << std::endl;
//    
//    // Print the results
//    std::cout << "Frequency: " << frequency << " Hz" << std::endl;
//   std::cout << "RMS value: " << rms << std::endl;
//
//    // Clean up
//    fftw_destroy_plan(plan);
//    fftw_cleanup_threads();
//    fftw_free(in);
//    fftw_free(out);
//
//    return 0;
//}
/*******************************************************************************
* Copyright 2015-2021 Intel Corporation.
*
* This software and the related documents are Intel copyrighted  materials,  and
* your use of  them is  governed by the  express license  under which  they were
* provided to you (License).  Unless the License provides otherwise, you may not
* use, modify, copy, publish, distribute,  disclose or transmit this software or
* the related documents without Intel's prior written permission.
*
* This software and the related documents  are provided as  is,  with no express
* or implied  warranties,  other  than those  that are  expressly stated  in the
* License.
*******************************************************************************/

//#include <stdio.h>
//#include "ipp.h"
//#include <chrono>
//#include <iostream>
//
///* Next two defines are created to simplify code reading and understanding */
//#define EXIT_MAIN exitLine:                                  /* Label for Exit */
//#define check_sts(st) if((st) != ippStsNoErr) goto exitLine; /* Go to Exit if Intel(R) Integrated Primitives (Intel(R) IPP) function returned status different from ippStsNoErr */
//
///* Results of ippMalloc() are not validated because Intel(R) IPP functions perform bad arguments check and will return an appropriate status  */
//
//int main()
//{
//    int len = 1024*16;
//    int i;
//    Ipp16u* pSrcDst = ippsMalloc_16u(len * sizeof(Ipp16u));
//    IppStatus status;
//    int scale = 0;  //without scaling
//    printf("\n\nSource vector\n");
//    for (i = 0; i < len; i++)
//    {
//        pSrcDst[i] = i;
//       // printf("%d; ", pSrcDst[i]);
//    }
//    
//    auto start = std::chrono::high_resolution_clock::now();
//
//    check_sts(status = ippsSqr_16u_ISfs(pSrcDst, len, scale));
//    /*
//    printf("\n\nResult\n");
//    for (i = 0; i < len; i++) printf("%d; ", pSrcDst[i]);
//    printf("\n\n");*/
//
//    EXIT_MAIN
//    auto end = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//    std::cout << "Time taken by function: " << duration.count() << "microsec" << std::endl;
//
//        ippsFree(pSrcDst);
//    printf("Exit status %d (%s)\n", (int)status, ippGetStatusString(status));
//    return (int)status;
//}
