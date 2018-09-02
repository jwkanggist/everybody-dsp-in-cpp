//
//  EverybodyDSPLib.cpp
//
//  Created by jwkangmacpro on 2017. 12. 7..
//  Copyright © 2017년 jwkangmacpro. All rights reserved.
//
#include "EverybodyDSPLib.h"


double sigmoid(const double a)
{
    return 1.0 / (1.0 + exp(-a) );
}


// function for raise-consine filter with roll-of factor
bool getRcwin(const unsigned int length, const double rolloff, const double Ts, double* rc_out)
{
    //postulation check
    
    if (length < 1)
    {
        LOGE("The length of raised cosine window should be positive.\n");
        return false;
    }
    else if (rolloff >1 || rolloff < 0)
    {
        LOGE("The roll-off factor should be less than or equal to '1' and more than or equal to '0'.\n");
        return false;
    }
    
    
    double time_stamp=-0.5*length*Ts;
    
    const double Trcosine = length*Ts;
    
    for (unsigned int i=0;i<length ;i++)
    {
        
        if (fabs(time_stamp) <= (1-rolloff)*Trcosine/(1+rolloff)/2.0 )
            rc_out[i]= 1.0;
        else if ( fabs(time_stamp) > Trcosine/2.0 )
            rc_out[i] =0.0;
        else
            rc_out[i] = 0.5*(1+cos((1.0+rolloff)*PI/rolloff/Trcosine*(fabs(time_stamp) - (1.0-rolloff)*Trcosine/(1.0+rolloff)/2.0)));
        
        time_stamp = time_stamp + Ts;
        
    }
    
    return true;
    
}

void LPFFIR(const float *FirCoeff, const unsigned int NumTaps, const double *Signal, double *FilteredSignal, const unsigned int NumSigPts)
{
    
    unsigned int j,  Top = 0;
    int k, n=0;
    double y, Reg[256];  // This assumes <= 256 taps.
    
    for (j = 0; j < NumTaps; j++)Reg[j] = 0.0;
    
    for (j = 0; j < NumSigPts + NumTaps; j++)
    {
        if (j < NumSigPts)
        {
            Reg[Top] = Signal[j];
            //LOGV("LPF: Signal[%d] = %lf", j, Signal[j]);
        }
        else
        {
            Reg[Top] = 0;
        }
        y = 0.0;
        n = 0;
        
        // The FirCoeff index increases while the Reg index decreases.
        for (k = Top; k >= 0; k--)
        {
            y += FirCoeff[n++] * Reg[k];
        }
        for (k = NumTaps - 1; k > (int)Top; k--)
        {
            y += FirCoeff[n++] * Reg[k];
        }
        FilteredSignal[j] = y;
        //LOGV("LPF: filteredSignal[%d] = %lf", j, FilteredSignal[j]);
        Top++;
        if (Top >= NumTaps)Top = 0;
    }
}


//----------------------------------------
// getPSDByFFT()
// operation: This function generate a power density spectrum [FFTSIZE/2+1 by 1] from input time domain signal inputTimeDomain [inputLen by 1] sampled by sampling rate Fs
//
// arguments:       1) double* inputtimeDomain: pointer of input time domain signal
//                  2) const unsigned int inputLen: the length of the input
//                  3) const unsigned int : sampling rate Fs
//                  4) const unsigned int FFTSIZE : FFTSIZE
//                  5) const unsigned int scalingConst : scaling constant for psd estimation
//                  6) double*  psd_out : pointer of the output psd
// final update @ 2017 12 11 Jaewook Kang
//----------------------------------
void getPSDByFFT(double* inputTimeDomain, const unsigned int inputLen, const unsigned int Fs, const unsigned int FFTSIZE, const double scalingConst, double* psd_out)
{
    EverybodyFFTLib* fftworker = new EverybodyFFTLib();
    
    // fft config
    fftworker->setRealFFT(true);
    fftworker->setFFTSIZE(FFTSIZE);
    fftworker->allocFFT(FFTSIZE);
    
    // Local memory allocation
    float_complex* inputFFTdomain = new float_complex[FFTSIZE];
    memset(inputFFTdomain,0,sizeof(float_complex)*FFTSIZE);
    
    float* inputTimeDomainFloat = new float[FFTSIZE];
    memset(inputTimeDomainFloat,0,sizeof(float)*FFTSIZE);
    
    double temp_ABSFFT=0.0;
    
    
    fftworker->fft(inputLen,inputTimeDomain,inputFFTdomain,inputTimeDomainFloat,inputFFTdomain);
    
    /*
     for (int ii=0;ii<5;ii++)
     {
     // cout<<"inputTimeDomain[i]="<<inputTimeDomain[ii]<<endl;
     cout<<"inputFFTdomain[i].r="<<inputFFTdomain[ii].r<<endl;
     cout<<"inputFFTdomain[i].i="<<inputFFTdomain[ii].i<<endl;
     }
     */
    // obtain power density spectrum
    
    for (unsigned int i=0 ; i < FFTSIZE/2+1 ;i++)
    {
        temp_ABSFFT=0.0;
        C_ABS(temp_ABSFFT,inputFFTdomain[i]);
        
        //cout<<"tempABS="<<temp_ABSFFT<<endl;
        
        if ( i==0 || i == FFTSIZE/2)
        {
            psd_out[i] = (double)((1.0 / scalingConst ) * pow(temp_ABSFFT,2));
            
        }
        else
        {
            psd_out[i] = (double)((2.0 / scalingConst) * pow(temp_ABSFFT,2));
            
        }
        
    }
    
    /*
     for (int i=0;i<5;i++)
     {
     cout<<"psd_out[i]="<<psd_out[i]<<endl;
     }
     */
    delete [] inputFFTdomain;
    delete [] inputTimeDomainFloat;
    
    inputFFTdomain = nullptr;
    inputTimeDomainFloat = nullptr;
    
    delete fftworker;
    
}


//------------------------------------------------
// iirSOSFilteringUnitSection()
// Operation:
// y[n] =  b01*x[n]   + b11*x[n-1] + b21*x[n-2]
//        -a11*y[n-1] - a21*y[n-2]
// where  SOS = [ b01 b11 b21 a01 a11 a21 ]
// Jwkang @ 2016 May
//------------------------------------------------
void iirSOSFilteringUnitSection(const double* IIRin, const unsigned int inputLen, const double* SOS, double* IIRout)
{
    double x[3]={0}; // input register  [ x[n] x[n-1] x[n-2] ]
    double y[3]={0}; // output register [ y[n] y[n-1] y[n-2] ]
    
    
    for (unsigned int i=0 ; i < inputLen; i++)
    {
        // shift register
        for (unsigned int n=2 ;n > 0 ;n--)
        {
            x[n]=x[n-1];
            y[n]=y[n-1];
        }
        
        // insert new input
        x[0]=IIRin[i];
        
        // calculation
        y[0] = SOS[0]*x[0] + SOS[1]*x[1] + SOS[2]*x[2] - SOS[4]*y[1] -SOS[5]*y[2];
        
        // generate new output
        IIRout[i]=y[0];
    }
    
}

//----------------------------------------
// Second order (biquadratic) IIR filtering
// iirSOSFilter()
// operation: This function filters the input data, IIRin, by sosIIRcoeff matrix.
//             the sosIIRcoeff matrix must be expressed using an [L x 6] second-order section (SOS)
//             matrix where L is the number of second-order sections.
//
//             The SOS matrix should have the following form:
//
//              SOS = [ b01 b11 b21 a01 a11 a21
//                      b02 b12 b22 a02 a12 a22
//                      ...
//                      b0L b1L b2L a0L a1L a2L ]
// arguments:       1) const float* IIRin : filter input
//                  2) const unsigned int inputLen: the length of the filter input
//                  3) double* sosIIRcoeff : IIR filter coefficient in the form of SOS
//                     sosIIRCoeff[i+iirOrder + j]

//                  4) const unsigned int sectionNum       : the number of sos sections
//                  5) double*                              : filter output
// Jaewook Kang @ 2016 May
//----------------------------------
void iirSOSFiltering(const double* IIRin, const unsigned int inputLen, const double* sosIIRCoeff, const unsigned int sectionNum, const unsigned int iirOrder, double* IIRout)
{
    double* reg = new double[inputLen];
    memcpy(reg, IIRin, inputLen*sizeof(double));
    
    for (unsigned int i=0; i < sectionNum ; i++)
    {
        iirSOSFilteringUnitSection(reg, inputLen, &sosIIRCoeff[i*iirOrder], IIRout);
        memcpy(reg,IIRout, inputLen*sizeof(double));
    }
    
    delete [] reg;
    reg = nullptr;
}


//----------------------------------------
// Filter the given input array using LPF
// filtering()
// operation: This function performs the filtering process with any given LPF profile pointed by inFirLpfProfilePtr
//
// arguments:       1) FIR_LPF_PROFILE *inFirLpfProfilePtr : The input pointer of FIR LPF profile struct
//                  2) double* inArray                     : The input data array
//                  3) const unsigned int nInArray         : The size of the input data array
//                  4) double* outArray                    : The output data array
//                  5) const bool deleteDelayFlag          : The flag on whether the fir filter delay is deleted or not
// Soonwon Ka @ 2017 Dec.
//----------------------------------
double* firFiltering(FIR_LPF_PROFILE *inFirLpfProfilePtr,
                  double* inArray, const unsigned int nInArray,
                  double* outArray,
                  const bool deleteDelayFlag) {
    double *ret;
    
    if(inFirLpfProfilePtr == nullptr) {
        // if the FIR LPF is not defined, do nothing and bypass the input array
        outArray = inArray;
        ret = outArray;
    }
    else {
        // lpf
        LPFFIR(inFirLpfProfilePtr->coeff, inFirLpfProfilePtr->numTaps,
               inArray, outArray, nInArray);
        
        // shifting for the right result
        if(deleteDelayFlag) {
            if(inFirLpfProfilePtr->numTaps%2 == 0) {
                ret = (double *)outArray+(inFirLpfProfilePtr->numTaps>>1);
            } else {
                ret = (double *)outArray+(inFirLpfProfilePtr->numTaps-1)/2;
            }
        } else {
            ret = (double *)outArray;
        }
    }
    // return the pointer of output array
    return ret;
}

bool downsampling(const size_t DOWN_SAMPLE_FACTOR,
                  const double* input, const size_t inputLen,
                  double* output, size_t outputLen)
{
    // postulation check
    if (DOWN_SAMPLE_FACTOR%2 != 0) return false;
    
    for(size_t i=0; i<inputLen;i++)
    {
        size_t idx = i / DOWN_SAMPLE_FACTOR;
        if(i%DOWN_SAMPLE_FACTOR== 0 && idx < outputLen) {
            
            output[idx] = input[i];
        }
    }
    
    return true;
}

bool downsampling(const size_t DOWN_SAMPLE_FACTOR, 
                  const double* input, const size_t inputLen, 
                  double* output, size_t outputLen,
                  const double scalingConst)
{
    // postulation check
    if (DOWN_SAMPLE_FACTOR%2 != 0) return false;
    
    for(size_t i=0; i<inputLen;i++)
    {
        size_t idx = i / DOWN_SAMPLE_FACTOR;
        if(i%DOWN_SAMPLE_FACTOR== 0 && idx < outputLen) {
            
            output[idx] = input[i]*scalingConst;
        }
    }
    
    return true;
}

