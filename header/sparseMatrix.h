//
//  sparseMatrix.h
//  sparsematrix
//
//  Created by jwkangmacpro on 2017. 5. 15..
//  Copyright © 2017년 jwkangmacpro. All rights reserved.
//  This codes is from jaewook kang's personal public GIT repo:
//  https://github.com/jwkanggist/sparseMatrixLib
//
// Designed by Jaewook Kang
// final update @ 2017 Oct

#ifndef sparseMatrix_h
#define sparseMatrix_h

#include <cstring> // for memset
#include <cstdio>
#include <cmath>
//#include <memory>
#include "log.h"



# ifdef __DSPLIBTEST__
#include <fstream> // not used in android
#include <iostream> // not used in android
#endif

const int INVAILD_VALUE = -999;

typedef struct SparseMatrixRepTable
{
    unsigned int rowIndex;
    unsigned int colIndex;
    double values;
    
} SparseMatrixRepTable;

class SparseMatrix{
private:

    SparseMatrixRepTable* mMtxTable = nullptr;
    
    const int mTableSize = 3;
    
    unsigned int mMatrixRowSize;
    unsigned int mMatrixColSize;
    unsigned int mNumOfNonZeros;
    
    void allocate(const unsigned  int numOfNonzeros_);
    void deAllocate(void);
    bool sortElement();

# ifdef __DSPLIBTEST__
    void showElemInMatrix(const int nonzeroIndex);
    void showElemInMatrixStructForm(const int nonzeroIndex);
#endif
    
    
public:
    SparseMatrix();
    SparseMatrix(const unsigned int rowSize_, const unsigned int colSize_);
    SparseMatrix(const unsigned int rowSize_, const unsigned int colSize_, const unsigned int nonZeroNum_);

    ~SparseMatrix();
    
    void setMatrix(SparseMatrixRepTable* table,const unsigned int rowSize_, const unsigned int colSize_, const unsigned int nonZeroNum_);
    
# ifdef __DSPLIBTEST__
    void readMatrix(const char* inputfilename);
    bool showMatrix(const char* str);
    bool showMatrixStructForm();
#endif
    
    unsigned int getRowSize() const { return mMatrixRowSize;}
    unsigned int getColSize() const { return mMatrixColSize;}
    unsigned int getNumofNonZeros() const { return  mNumOfNonZeros;}
    
    int getRowIndex(const int nonzeroIndex) const
    {
        if ( nonzeroIndex < mNumOfNonZeros)
            return mMtxTable[nonzeroIndex].rowIndex;
        else
            return -1;
    }
    
    
    int getColIndex(const int nonzeroIndex) const
    {
        if ( nonzeroIndex < mNumOfNonZeros)
            return mMtxTable[nonzeroIndex].colIndex;
        else
            return -1;
    }
    
    double getElemValue(const int nonzeroIndex) const
    {
        if (mMtxTable != nullptr)
        {
            if ( nonzeroIndex < mNumOfNonZeros)
                return mMtxTable[nonzeroIndex].values;
            else
                return 0.0;
        }
        else
        {
            LOGE("[SparseMatrix] mMtxTable is not allocated!");
            return INVAILD_VALUE;
        }
    }
    
    double getElemValue(const unsigned int rowIndex, const unsigned int colIndex);
    
    SparseMatrixRepTable* getSparseMatrixRepTable() const { return mMtxTable;}

    bool transpose(SparseMatrix* X);
    bool transpose();
    bool isAllocated() const {
        if (mMtxTable == nullptr)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    
    

};



#endif /* sparseMatrix_h */
