//
//  sparseMatrix.cpp
//  sparsematrix
//
//  Created by jwkangmacpro on 2017. 5. 15..
//  Copyright © 2017년 jwkangmacpro. All rights reserved.
//
//  This codes is from jaewook kang's personal public GIT repo:
//  https://github.com/jwkanggist/sparseMatrixLib
//
//

#include "sparseMatrix.h"


SparseMatrix::SparseMatrix()
{
    mMatrixRowSize = 0;
    mMatrixColSize = 0;
    mNumOfNonZeros = 0;
# ifdef __DSPLIBTEST__
    std::cout<<"------------------------------"<<std::endl;
    std::cout<<"[sparseMatrix] Matrix init"<<std::endl;
    std::cout<<"[sparseMatrix] mMatrixRowSize = "<<mMatrixRowSize<<std::endl;
    std::cout<<"[sparseMatirx] mMatrixColSize = "<<mMatrixColSize<<std::endl;
#endif
}


SparseMatrix::SparseMatrix(const  unsigned int rowSize_, const  unsigned int colSize_)
{
    mMatrixRowSize = rowSize_;
    mMatrixColSize = colSize_;
    mNumOfNonZeros = 0;
# ifdef __DSPLIBTEST__
//    std::cout<<"------------------------------"<<std::endl;
//    std::cout<<"[sparseMatrix] Matrix init"<<std::endl;
//    std::cout<<"[sparseMatrix] mMatrixRowSize = "<<mMatrixRowSize<<std::endl;
//    std::cout<<"[sparseMatirx] mMatrixColSize = "<<mMatrixColSize<<std::endl;
#endif

}

SparseMatrix::SparseMatrix(const  unsigned int rowSize_, const  unsigned int colSize_, const  unsigned int nonZeroNum_)
{
    mMatrixRowSize = rowSize_;
    mMatrixColSize = colSize_;
    mNumOfNonZeros = nonZeroNum_;
# ifdef __DSPLIBTEST__
    std::cout<<"------------------------------"<<std::endl;
    std::cout<<"[sparseMatrix] Matrix init"<<std::endl;
    std::cout<<"[sparseMatrix] mMatrixRowSize = "<<mMatrixRowSize<<std::endl;
    std::cout<<"[sparseMatirx] mMatrixColSize = "<<mMatrixColSize<<std::endl;
    std::cout<<"[sparseMatirx] mNumOfNonZeros = "<<mNumOfNonZeros<<std::endl;
#endif
}

SparseMatrix::~SparseMatrix()
{
    // memory deallocation
    if ( mMtxTable != nullptr)
    {
        deAllocate();
    }
}


void SparseMatrix::allocate(const  unsigned int numOfNonzeros_)
{
    // memory allocation and init
    // -- allocating memory by considering the maximum case of the nonzeros in the sparse matrix.

    mMtxTable = new SparseMatrixRepTable[numOfNonzeros_];
    memset(mMtxTable,0.0,sizeof(SparseMatrixRepTable)*numOfNonzeros_);
    
   // std::cout<<"[sparseMatrix] Allocate [unsigned int: 2 x "<<numOfNonzeros_<<"] memory for matrix storage"<<std::endl;
   // std::cout<<"[sparseMatrix] Allocate [double: 1 x "<<numOfNonzeros_<<"] memory for matrix storage"<<std::endl;
}

void SparseMatrix::deAllocate(void)
{
    if (mMtxTable !=nullptr)
    {
        delete [] mMtxTable;
        mMtxTable = nullptr;
    }
    
//    std::cout<<"------------------------------"<<std::endl;
  //  std::cout<<"[sparseMatrix] deAllocate memory"<<std::endl;
    //std::cout<<"--- mMtxTable="<<mMtxTable<<std::endl;

}

# ifdef __DSPLIBTEST__
void SparseMatrix::readMatrix(const char* inputfilename)
{
    SparseMatrixRepTable* tempSparseMatrixRepTable = nullptr;
    
    tempSparseMatrixRepTable = new SparseMatrixRepTable[mMatrixColSize*mMatrixRowSize];
    memset(tempSparseMatrixRepTable, 0.0, sizeof(SparseMatrixRepTable)*mMatrixColSize*mMatrixRowSize);
    
    std::ifstream inputfile;
    double tempReadValue = 0;
    
    inputfile.open(inputfilename);
    unsigned int nonZerosCnt = 0;

    std::cout<<"[sparseMatrix] File read from "<<inputfilename<<std::endl;
    std::cout<<"------------------------------"<<std::endl;
    std::cout<<"[sparseMatrix] Storing the sparse matrix via the triple representation method"<<std::endl;
    std::cout<<"Imported Matrix ="<<std::endl;
    for (int rowIndex = 0 ; rowIndex < mMatrixRowSize ; rowIndex++)
    {
        //std::cout<<"------- "<<rowIndex<<"-th row"<<" ----------"<<std::endl;
        for (int colIndex = 0 ; colIndex < mMatrixColSize ; colIndex++)
        {
            inputfile >> tempReadValue;
            //std::cout << "[sparseMatrix] tempReadValue = " <<  tempReadValue << std::endl;
            std::cout<<tempReadValue<<",";
            if (tempReadValue != 0)
            {
                tempSparseMatrixRepTable[nonZerosCnt].rowIndex = rowIndex;
                tempSparseMatrixRepTable[nonZerosCnt].colIndex = colIndex;
                tempSparseMatrixRepTable[nonZerosCnt].values = tempReadValue;
                nonZerosCnt++;
            }
        }
        std::cout<<std::endl;
    }
    inputfile.close();
    
    mNumOfNonZeros = nonZerosCnt;
    // memcpy and free the temporal memory
    allocate(mNumOfNonZeros);
    
    memcpy(mMtxTable,tempSparseMatrixRepTable,sizeof(SparseMatrixRepTable)*mNumOfNonZeros);
    
    delete [] tempSparseMatrixRepTable;
    tempSparseMatrixRepTable =nullptr;
    
    
    std::cout<<"------------------------------"<<std::endl;
    std::cout<<"--- Total matrix size: "<<mMatrixRowSize<<" X "<< mMatrixColSize<<" = "<<mMatrixRowSize*mMatrixColSize<<std::endl;
    std::cout<<"--- Total number of nonzeros = "<<nonZerosCnt<<std::endl;
    std::cout<<"--- Matrix sparsity = "<< (float)nonZerosCnt / (float)(mMatrixRowSize*mMatrixColSize) * 100.0 <<" %"<<std::endl;
    std::cout<<"--- Total number of memory slots for triple representation = 3 X "<<nonZerosCnt<<" = "<<nonZerosCnt*mTableSize<<std::endl;
    
}



bool SparseMatrix::showMatrix(const char* str)
{
    if (mMtxTable !=nullptr)
    {
    //    std::cout<<"------------------------------"<<std::endl;
        std::cout<<"[sparseMatrix] The triple representation of the sparse matrix"<<std::endl;
        std::cout<<"mMatrixRowSize="<<mMatrixRowSize<<",mMatrixColSize="<<mMatrixColSize<<std::endl;
        std::cout<<"mNumOfNonZeros="<<mNumOfNonZeros<<std::endl;

        std::cout<<"(rowIndex,colIndex), nonzerovalue"<<std::endl;
        
        std::cout<<str<<"="<<std::endl;
        for (int elemIndex = 0 ; elemIndex < mNumOfNonZeros ; elemIndex ++)
        {
            showElemInMatrix(elemIndex);
        }
        return true;
    }
    else
    {
        LOGD("[SparseMatrix] mMtxTable is not allocated!");
        return false;
    }
}

bool SparseMatrix::showMatrixStructForm()
{
    if (mMtxTable !=nullptr)
    {
        //    std::cout<<"------------------------------"<<std::endl;
        std::cout<<"[sparseMatrix] The triple representation of the sparse matrix"<<std::endl;
        std::cout<<"(rowIndex,colIndex), nonzerovalue"<<std::endl;
        
        for (int elemIndex = 0 ; elemIndex < mNumOfNonZeros ; elemIndex ++)
        {
            showElemInMatrixStructForm(elemIndex);
        }
        return true;
    }
    else
    {
        LOGD("[SparseMatrix] mMtxTable is not allocated!");
        return false;
    }
}

void SparseMatrix::showElemInMatrix(const int nonzeroIndex)
{
    std::cout<<"("<<mMtxTable[nonzeroIndex].rowIndex<<","<<mMtxTable[nonzeroIndex].colIndex<<"),"<<" "<<mMtxTable[nonzeroIndex].values<<std::endl;
}

void SparseMatrix::showElemInMatrixStructForm(const int nonzeroIndex)
{
    std::cout<<"{"<<mMtxTable[nonzeroIndex].rowIndex<<","<<mMtxTable[nonzeroIndex].colIndex<<","<<mMtxTable[nonzeroIndex].values<<"},"<<std::endl;
}




#endif

double SparseMatrix::getElemValue(const unsigned int rowIndex, const unsigned int colIndex)
{
    if (mMtxTable != nullptr)
    {
        int elemIndex = -1;


        for (int i = 0 ; i < mNumOfNonZeros ; i++)
        {
            if (mMtxTable[i].colIndex == colIndex && mMtxTable[i].rowIndex == rowIndex)
            {
                elemIndex = i;
                break;
            }
        }
        
        if (elemIndex < 0)
        {
            return 0.0;
        }
        else
        {
            // showElemInMatrix(elemIndex);
            return mMtxTable[elemIndex].values;
        }
    }
    else
    {
        LOGE("[SparseMatrix] mMtxTable is not allocated!");
        return INVAILD_VALUE;
    }
    
}






void SparseMatrix::setMatrix(SparseMatrixRepTable* table,const unsigned int rowSize_, const unsigned int colSize_, const unsigned int nonZeroNum_)
{
    mMatrixRowSize = rowSize_;
    mMatrixColSize = colSize_;
    mNumOfNonZeros = nonZeroNum_;
    
    allocate(mNumOfNonZeros);
    memcpy(mMtxTable,table,sizeof(SparseMatrixRepTable)*mNumOfNonZeros);
    
}

bool SparseMatrix::transpose(SparseMatrix* X)
{
    if (X->isAllocated())
    {
        mMatrixColSize = X->getRowSize();
        mMatrixRowSize = X->getColSize();
        mNumOfNonZeros = X->getNumofNonZeros();
        
        allocate(mNumOfNonZeros);
        for (unsigned int n=0 ; n < mNumOfNonZeros; n++)
        {
            mMtxTable[n].rowIndex = X->getColIndex(n);
            mMtxTable[n].colIndex = X->getRowIndex(n);
            mMtxTable[n].values   = X->getElemValue(n);
        }
        sortElement();
        return true;
    }
    else
    {
        LOGE("[SparseMatrix] X->mMtxTable is not allocated!");
        return false;
    }
}

bool SparseMatrix::transpose()
{
    if ( mMtxTable !=nullptr)
    {
        unsigned int tempColSize = getColSize();
        
        mMatrixColSize = getRowSize();
        mMatrixRowSize = tempColSize;
        
        for (unsigned int n=0 ; n < mNumOfNonZeros; n++)
        {
            unsigned int tempColIndex = getColIndex(n);
            mMtxTable[n].colIndex = getRowIndex(n);
            mMtxTable[n].rowIndex = tempColIndex;
        }
        sortElement();
        return true;
    }
    else
    {
        LOGE("[SparseMatrix] mMtxTable is not allocated!");
        return false;
    }
}

bool SparseMatrix::sortElement()
{
    if (mMtxTable != nullptr)
    {
        unsigned int m = 0;
        for (unsigned int n = 1; n < mNumOfNonZeros; n++)
        {
            unsigned int tempRowIndex = getRowIndex(n);
            SparseMatrixRepTable tempTable = mMtxTable[n];
            m = n;
            while ( getRowIndex(m-1) > tempRowIndex && m > 0)
            {
                mMtxTable[m] = mMtxTable[m-1];
                m--;
            }
            mMtxTable[m] = tempTable;
        }
        return true;
    }
    else
    {
        LOGE("[SparseMatrix] mMtxTable is not allocated!");
        return false;
    }
}
