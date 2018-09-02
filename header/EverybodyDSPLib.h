//
//  EverybodyDSPLib.h
//
//  Created by jwkangmacpro on 2017. 12. 7..
//  Copyright © 2017년 jwkangmacpro. All rights reserved.
//

#ifndef EverybodyDSPLib_h
#define EverybodyDSPLib_h
#include "log.h"
#include <cmath>
#include "EverybodyMathLib.h"
#include "EverybodyFFTLib.h"
#include "JacobiSvd.h"

struct FIR_LPF_PROFILE
{
    float *coeff=nullptr;
    unsigned int numTaps;
    FIR_LPF_PROFILE();
    FIR_LPF_PROFILE(float *inCoeff, const unsigned int inNumTaps) {
        coeff = inCoeff;
        numTaps = inNumTaps;
    }
};


double sigmoid(const double a);

// Generation of raised cosine window //
bool getRcwin(const unsigned int length, const double rolloff, const double Ts, double* rc_out);


// FIR filtering for LPF
// HPF or BPF is also compatible if an appropriate coefficient is prepared.
void LPFFIR(const float *FirCoeff, const unsigned int NumTaps,
            const double *Signal, double *FilteredSignal, const unsigned int NumSigPts);

// IIR filtering using second order section structure
void iirSOSFilteringUnitSection(const double* IIRin, const unsigned int inputLen, const double* SOS, double* IIRout);
void iirSOSFiltering(const double* IIRin, const unsigned int inputLen, const double* sosIIRCoeff, const unsigned int sectionNum, const unsigned int iirOrder, double* IIRout);

void getPSDByFFT(double* inputTimeDomain, const unsigned int inputLen, const unsigned int Fs, const unsigned int FFTSIZE, const double scalingConst, double* psd_out);


double* firFiltering(FIR_LPF_PROFILE *inFirLpfProfilePtr, double* inArray, const unsigned int nInArray, double* outArray, const bool deleteDelayFlag);

bool downsampling(const size_t DOWN_SAMPLE_FACTOR, const double* input, const size_t inputLen, double* output, size_t outputLen);
bool downsampling(const size_t DOWN_SAMPLE_FACTOR, const double* input, const size_t inputLen, double* output, size_t outputLen, const double scalingConst);


#endif /* EverybodyDSPLib_h */
