//#include <iostream>
//#include <vector>
//#include <complex>
//#include <fftw3.h>
//#include <omp.h>
//#include <experimental/vector>
//#include <chrono>
//#include <atomic>
//#include <algorithm>
//#include <ipp.h>
//
//const double PI = 3.1415926535897932384;
//#define EXIT_MAIN exitLine:                                  /* Label for Exit */
//#define check_sts(st) if((st) != ippStsNoErr) goto exitLine; 
//
//int main()
//{
//	const int n = 1024 * 16;
//	Ipp64f* pSRC = ippsMalloc_64f(n * sizeof(Ipp64f));
//	Ipp64f* pDST = ippsMalloc_64f(n * sizeof(Ipp64f));
//	Ipp64f* mag  = ippsMalloc_64f(n * sizeof(Ipp64f));
//	Ipp8u* pMemInit = NULL, * pBuffer = NULL, * pSpecMem = NULL; /* Pointer to the work buffers */
//	int sizeSpec = 0, sizeInit = 0, sizeBuf = 0;               /* size of FFT pSpec structure, Init and work buffers */
//	IppStatus status = ippStsNoErr;
//	const double fs = 100000;
//	//    // Sampling interval
//	    const double T = 1.0 / fs;
//	//    // Wave frequency
//	    const double f = 100;
//		double t;
//
//		
//		check_sts(status = ippsFFTGetSize_R_64f(16, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, &sizeSpec, &sizeInit, &sizeBuf))
//		pSpecMem = (Ipp8u*)ippMalloc(sizeSpec);
//		pBuffer = (Ipp8u*)ippMalloc(sizeBuf);
//		pMemInit = (Ipp8u*)ippMalloc(sizeInit);
//		check_sts(status = ippsFFTInit_R_64f(&pSpec, 16, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, pSpecMem, pMemInit))
//		
//			for (int i = 0; i < n; i++)
//			{
//				t = i * T;
//				pSRC[i] = sin(2 * PI * f * t);
//
//			}
//
//		check_sts(status = ippsFFTFwd_RToCCS_64f(pSrc, pDst, pSpec, pBuffer))
//
//		check_sts(status = ippsMagnitude_64f((Ipp64fc*)pDst, pMagn, n))
//
//EXIT_MAIN
//		ippFree(pMemInit);
//		ippFree(pSpec);
//		ippFree(pBuffer);
//		return (int)status;
//
//}