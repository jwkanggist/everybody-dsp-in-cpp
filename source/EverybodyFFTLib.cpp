//
//  EverybodyFFTLib.cpp
//
//  Created by jwkangmacpro on 2017. 10. 13..
//  Copyright © 2017년 jwkangmacpro. All rights reserved.
//


#include "EverybodyFFTLib.h"

EverybodyFFTLib::EverybodyFFTLib(){}


EverybodyFFTLib::EverybodyFFTLib(ConstParams* inConstParams, const unsigned int fftsize)
{
    //LOGD("> call EverybodyFFTLib()");

    setRealFFT(inConstParams->isRealFFT());
    setFFTSIZE(fftsize);
    allocFFT(fftsize);
    
}

EverybodyFFTLib::~EverybodyFFTLib()
{
    //LOGD("> call ~EverybodyFFTLib()");

    deallocFFT();
}

void EverybodyFFTLib::absForcpx(const unsigned length, float_complex* cpxIn, float* absout)
{
    for (unsigned int n = 0 ; n < length ; n++)
    {
        C_ABS(absout[n],cpxIn[n]);
    }
    
}

void EverybodyFFTLib::phaseForcpx(const unsigned length, float_complex* cpxIn, float* phaseout)
{
    for (unsigned int n = 0 ; n < length ; n++)
    {
        phaseout[n] = atan2f(cpxIn[n].i,cpxIn[n].r);
    }
    
}

void EverybodyFFTLib::real2cpx(const unsigned length, const double *real, float_complex *cpx)
{
    for(unsigned int i = 0; i < length; i++)
    {
        cpx[i].r = (double)real[i];
    }
}
 
void EverybodyFFTLib::imag2cpx(const unsigned length, const double *imag, float_complex *cpx)
{
    for(unsigned int i = 0; i < length; i++)
    {
        cpx[i].i = (double)imag[i];
    }
}

void EverybodyFFTLib::double2float(const unsigned length, const double *input, float* output)
{
    for( unsigned int i=0; i< length; i ++ )
    {
        output[i] = (float)input[i];
    }
}

void EverybodyFFTLib::fft(const unsigned int minArrayLen,double *lfIn, float_complex* cpxIn, float* fIn, float_complex* cpxOut )
{
     if(mRealFFTOn == false)
     {
         memset(cpxIn, 0, mFftSize * sizeof(float_complex));
         real2cpx(minArrayLen, lfIn, cpxIn); // type conversion
         kiss_fft(mForwardFftCfg, cpxIn, cpxOut); // fft
     }
     else
     {
         memset(fIn, 0, mFftSize * sizeof(float));
         double2float(minArrayLen, lfIn, fIn);// type conversion from double to float (2017 Feb, swka)
         kiss_fftr(mForwardFftrCfg, fIn, cpxOut); // fft
         for (unsigned int i = 0; i < (mFftSize >> 1)-1; i++)
         {
             cpxOut[(mFftSize>>1)+1+i].r = cpxOut[(mFftSize>>1)-1-i].r;
             cpxOut[(mFftSize>>1)+1+i].i = (-1)*cpxOut[(mFftSize>>1)-1-i].i;
         }
     }
}


void EverybodyFFTLib::ifft(float_complex* cpxIn, float_complex* cpxOut )
{
    kiss_fft(mInverseFftCfg, cpxIn, cpxOut);

}

void EverybodyFFTLib::hilbertMul(float_complex* cpxOut, float_complex* in1, float_complex* in2)
{
    for (unsigned int i = (mFftSize >> 1); i < mFftSize; i++)
    {
        C_MUL(cpxOut[i], in1[i], in2[i]);
 //        printf("   output[%u]: %f,  in1[%u]: %f + j%f,  in2[%u]: %f + j%f\n",i, output[i], i, in1[i].r, in1[i].i, i, in2[i].r, in2[i].i);
    }
}


void EverybodyFFTLib::magPhaseCombine2Cpx(const unsigned int length, float* magArray, float* phaseArray, float_complex* cpxArrayOut)
{
    for (unsigned int n = 0 ; n < length ; n++)
    {
        float currphase = phaseArray[n];
        float currmag   = magArray[n];
        
        cpxArrayOut[n].r = currmag * cos(currphase);
        cpxArrayOut[n].i = currmag * sin(currphase);
//        LOGD("JammingCanceller] cpxArrayOut[%d].r=%e",n,cpxArrayOut[n].r);
//        LOGD("JammingCanceller] cpxArrayOut[%d].i=%e",n,cpxArrayOut[n].i);
//        LOGD("JammingCanceller] currphase[%d].r=%e",n,currphase);

    }
    
}


void EverybodyFFTLib::allocFFT(const unsigned int fftsize)
{
    mFftSize = fftsize;
    deallocFFT();

    if( mRealFFTOn == false)
    {
        mForwardFftCfg = kiss_fft_alloc(mFftSize, 0, NULL, NULL);
    }
    else
    {
        mForwardFftrCfg = kiss_fftr_alloc(mFftSize, 0, NULL, NULL);
    }
    mInverseFftCfg = kiss_fft_alloc(mFftSize, 1, NULL, NULL);
}

void EverybodyFFTLib::deallocFFT()
{
    if( mRealFFTOn == false)
    {
        if (mForwardFftCfg != NULL)
        {
            free(mForwardFftCfg);
            mForwardFftCfg = NULL;
        }

    }
    else
    {
        if (mForwardFftrCfg != NULL)
        {
            free(mForwardFftrCfg);
            mForwardFftrCfg = NULL;
        }

    }
    if (mInverseFftCfg != NULL)
    {
        free(mInverseFftCfg);
        mInverseFftCfg = NULL;
    }

}
