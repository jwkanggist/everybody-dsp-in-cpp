//
//  EverybodyFFTLib.h
//
//  Created by jwkangmacpro on 2017. 10. 13..
//  Copyright © 2017년 jwkangmacpro. All rights reserved.
//
//

#ifndef EverybodyFFTLib_h
#define EverybodyFFTLib_h



#include "_kiss_fft_guts.h"// we are using kiss_fft library (2017 Oct )


#include "tempConstParamsForTest.h"


#include <cmath>

// Type of the data structure for fft I/O
typedef kiss_fft_cpx float_complex;

class EverybodyFFTLib{
#ifndef __CORE__
protected:
#else
#ifdef UNIT_TEST
public:
#else
protected:
#endif
#endif
    bool mRealFFTOn;
    
    unsigned int mFftSize; // FFT size
    
    kiss_fft_cfg mForwardFftCfg = NULL; // kiss FFT instance
    kiss_fft_cfg mInverseFftCfg = NULL; // kiss IFFT instance
    kiss_fftr_cfg mForwardFftrCfg = NULL; // kiss real FFT instance (swka, Feb. 2017)
    
public:
    EverybodyFFTLib();
    EverybodyFFTLib(ConstParams* inConstParams,unsigned int fftsize);
    
    ~EverybodyFFTLib();
    
    //------------------------------------
    // set method
    void setRealFFT(const bool IN_REAL_FFT) {mRealFFTOn = IN_REAL_FFT;}
    void setFFTSIZE(const unsigned int FftSize){mFftSize = FftSize;}
    
    void allocFFT(const unsigned int fftsize);
    void deallocFFT();
    
    // working for
    // 1) complex-valued fft and (setRealFFT(false))
    // 2) real-valued fft both   (setRealFFT(true))
    void fft (const unsigned int minArrayLen,double *lfIn, float_complex* cpxIn, float* fIn, float_complex* output );
    
    // working only by complex-valued ifft
    void ifft(float_complex* cpxIn, float_complex* cpxOut );
    
    // some auxiliay functions for complex value manipulation
    void real2cpx(const unsigned length, const double *real, float_complex *cpx);
    void imag2cpx(const unsigned length, const double *imag, float_complex *cpx);
    void absForcpx  (const unsigned length, float_complex* cpxIn, float* absout); // return magnitude of cpxIn
    void phaseForcpx(const unsigned length, float_complex* cpxIn, float* phaseout);// return phase in radian
    void magPhaseCombine2Cpx(const unsigned int length, float* magArray, float* phaseArray, float_complex* cpxArrayOut);

    
    // for double to float array type conversion
    void double2float(const unsigned length, const double *input, float* output);
    void hilbertMul(float_complex* output, float_complex* in1, float_complex* in2);

    
    bool isRealFFT() const {return mRealFFTOn;}

    
};
#endif /* EverybodyFFTLib_h */
