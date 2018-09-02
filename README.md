
# Everybody DSP Library #
Author: Jaewook Kang

## About
Welcome to **Everybvody DSP in C++** !

The source codes in this repo are mainly have focused on efficient implementation of  linear algebra and DSP operation  for mobile platforms such as Android and iOS. 

This library include the below modules:
- `EverybodyDSPLib`:  which is a wrapper module containing the below two modules, additionally including several DSP methods.
- `EverybodyMathLib`: which includes C++ classes of Matrix and Vector, and some related functions for linear algebra manipulations. 
- `EverybodyFFTLib`: which is a wrapper module for an external FFT package (currently kissFFT), including some methods to manipulate complex number system. 


## Dependencies ##
- [kiss FFT](https://sourceforge.net/projects/kissfft/)
- [Rosetta Code's QRdecomp](https://rosettacode.org/wiki/QR_decomposition)
- [JwKang's Sparse Matrix repo](https://github.com/jwkanggist/sparseMatrixLib)

## How to use
You can submodule this repo in your other repository as following
```bash
git clone 
git init
```

For use, just include the header file of `EverybodytDSPLib.h` in your C++ codes.
```cpp
#include "EverybodyDSPLib.h"

```

## Features ##
#### v1.0
  * FFT operation modules and periphery functions (kiss FFT)
  * Estimation of power spectral density using FFT
  * Type-1 FIR filtering using linear convolution 
  * IIR filtering using second order section (SOS) structure
  * 1D downsampling  
  * Generation of raised cosine window
  * Matrix and Vector data structures
  * Some linear algebra operations
  * QR decomposition using householder transform
  * Jacobian singular value decomposition

## Test Platforms ##
  * Android <= 8.0 (Android Oreo)
  * iOS <= 11.0 (iOS11)
 
## JNI Interface for Android
- TBU

## Benchmark
- TBU

## Contributor Info
* Jaewook Kang (jwkang10@gmail.com)



