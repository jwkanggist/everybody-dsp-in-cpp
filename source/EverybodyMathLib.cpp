//
//  EverybodyMathLib.cpp
//
//  Created by jwkangmacpro on 2017. 4. 18..
//  Copyright © 2017년 jwkangmacpro. All rights reserved.
//
//  reference :  https://rosettacode.org/wiki/QR_decomposition,
//  some revision are done by jwkang 2017 Apr.

#include "EverybodyMathLib.h"

void Vector::powElem(const double exponent)
{
    for (unsigned int i = 0 ; i <size ; i++)
    {
        (*this)(i) = pow(data[i],exponent);
    }
}



double Vector::sumAllElem()
{
    double sum = 0.0;
    
    for (unsigned int i = 0 ; i < size ; i++ )
    {
        sum += data[i];
    }
    return sum;
}



double Vector::powSumAllElem(const double exponent)
{
    double sqrsum = 0.0;
    
    for (unsigned int i = 0 ; i < size ; i++ )
    {
        sqrsum += pow(data[i],exponent);
    }
    return sqrsum;
}



double Vector::absPowSumAllElem(const double exponent)
{
    double abspowsum = 0.0;
    for ( unsigned int i =0 ; i <size ; i++)
    {
        abspowsum += pow(fabs( data[i] ),exponent );
    }
    return abspowsum;
}




void Vector::absElem()
{
    for (unsigned int i = 0 ; i < size ; i++)
    {
        (*this)(i) = fabs(data[i]);
    }
}

void Vector::absSqrtElem()
{
    for (unsigned int i = 0 ;  i < size ; i++)
    {

        (*this)(i) = sqrt(fabs(data[i]));
    }
}




void Vector::absPowElem(const double exponent)
{
    for (unsigned int i = 0 ; i < size ; i++)
    {
        (*this)(i) = pow(fabs(data[i]),exponent);
    }
}




void Matrix::powElem(const double exponent)
{
    for (unsigned int i = 0 ;  i < m ; i++)
    {
        for (unsigned int j = 0 ; j < n ; j++)
        {
            (*this)(i,j) = pow(data[i+m*j],exponent);
        }
    }
}



void Matrix::absSqrtElem()
{
    for (unsigned int i = 0 ;  i < m ; i++)
    {
        for (unsigned int j = 0 ; j < n ; j++)
        {
            (*this)(i,j) = sqrt(fabs(data[i+m*j]));
        }
    }
}



double Matrix::powSumAllElem(const double exponent)
{
    double sqrsum = 0.0;
    for (unsigned int i = 0 ;  i < m ; i++)
    {
        for (unsigned int j = 0 ; j < n ; j++)
        {
            sqrsum +=  pow(data[i+m*j],exponent);
        }
    }
    
    return sqrsum;
}




double Matrix::sumAllElem()
{
    double sum = 0.0;
    for (unsigned int i = 0 ;  i < m ; i++)
    {
        for (unsigned int j = 0 ; j < n ; j++)
        {
            sum +=  data[i+m*j];
        }
    }
    
    return sum;
}




void Matrix::absElem()
{
    for (unsigned int i = 0 ;  i < m ; i++)
    {
        for (unsigned int j = 0 ; j < n ; j++)
        {
            (*this)(i,j) = fabs(data[i+m*j]);
        }
    }
}




double Matrix::absSumAllElem()
{
    double abssum = 0.0;
    for (unsigned int i = 0 ;  i < m ; i++)
    {
        for (unsigned int j = 0 ; j < n ; j++)
        {
            abssum +=  fabs(data[i+m*j]);
        }
    }
    
    return abssum;
}




// c = a + b * s
void vmadd(const Vector& a, const Vector& b, double s, Vector& c)
{
    if (c.size != a.size or c.size != b.size) {
        LOGE("[EverybodyMathLib] [vmadd]: vector sizes don't match");
        return;
    }
    
    for (int i = 0; i < c.size; i++)
        c(i) = a(i) + s * b(i);
}




// c = c1 * a + c2 * b
void vmadd(const Vector& a, const Vector& b, const double c1, const double c2, Vector& vecout)
{
    if (a.size != b.size) {
        LOGE("[EverybodyMathLib] [vmadd]: vector sizes don't match");
        return;
    }
    else
    {
        vecout.allocate(a.size);
        for (int i = 0; i < vecout.size; i++)
            vecout(i) = c1 * a(i) + c2 * b(i);
    }
    
}




// mat = I - 2*v*v^T
// !!! m is allocated here !!!
void compute_householder_factor(Matrix& mat, const Vector& v)
{
    
    int n = v.size;
    mat.allocate(n,n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            mat(i,j) = -2 *  v(i) * v(j);
    for (int i = 0; i < n; i++)
        mat(i,i) += 1;
}




// take c-th column of a matrix, put results in Vector v
void Matrix::extract_column(Vector& v_extracted, unsigned int c)
{
    /*
    if (m != v.size) {
        LOGE("[EverybodyMathLib] [Matrix::extract_column]: Matrix and Vector sizes don't match");
        return;
    }
    */
    v_extracted.allocate(m);
    for (int i = 0; i < m; i++)
        v_extracted(i) = (*this)(i,c);
}



// insert a vector into the c-th column of this matrix
void Matrix:: insert_column(Vector& v_inserted, unsigned int c)
{
    if (this->data != nullptr)
    {
        if ( (this->m) == v_inserted.size)
        {
            for (int i =0 ; i < m; i++)
            {
                (*this)(i,c) = v_inserted(i);
            }
        }
        else
        {
            LOGD("[Matrix] teh column size mismatch!");
        }
    }
    else
    {
        LOGD("[Matrix] This matrix is not allocated.");
    }
    
}




# ifdef __DSPLIBTEST__
void Matrix::readMatrix(const char* inputfilename)
{
    std::ifstream inputfile;
    double tempReadValue = 0.0;
    
    inputfile.open(inputfilename);
    std::cout<<"[Matrix] File read from "<<inputfilename<<std::endl;
    std::cout<<"------------------------------"<<std::endl;
    std::cout<<"Imported Matrix ="<<std::endl;
    
    for (unsigned int rowIndex = 0 ; rowIndex < m ; rowIndex++)
    {
        //std::cout<<"------- "<<rowIndex<<"-th row"<<" ----------"<<std::endl;
        for (unsigned int colIndex = 0 ; colIndex < n ; colIndex++)
        {
            inputfile >> tempReadValue;
            //std::cout << "[sparseMatrix] tempReadValue = " <<  tempReadValue << std::endl;
            data[rowIndex+colIndex*m] = tempReadValue;
            std::cout<<tempReadValue<<",";
        }
        std::cout<<std::endl;
    }
    inputfile.close();
    
    std::cout<<"------------------------------"<<std::endl;
    std::cout<<"--- Total matrix size: "<<m<<" X "<< n<<" = "<<m * n <<std::endl;
    
}



void Vector::readVector(const char* inputfilename)
{
    std::ifstream inputfile;
    double tempReadValue = 0.0;
    
    inputfile.open(inputfilename);
    std::cout<<"[Vector] File read from "<<inputfilename<<std::endl;
    std::cout<<"------------------------------"<<std::endl;
    std::cout<<"Imported Vector ="<<std::endl;
    
    for (unsigned int index = 0 ; index < size ; index++)
    {
        inputfile >> tempReadValue;
        data[index] = tempReadValue;
        std::cout<<tempReadValue<<",";
    }
    std::cout<<std::endl;

    inputfile.close();
    
    std::cout<<"------------------------------"<<std::endl;
    std::cout<<"--- Total vector size: "<<size<<" X 1"<<std::endl;
    
}

void Matrix::showMatrix(const char* str)
{
    std::cout<<"-------------------------------"<<std::endl;
    std::cout<<"mMatrixRowSize = "<<this->m<<std::endl;
    std::cout<<"mMatrixColSize = "<<this->n<<std::endl;
    
    std::cout<<str<<"="<<std::endl;
    for(unsigned int i = 0; i < this->m; i++) {
        std::cout<<"{";
        for (unsigned int j = 0; j < this->n; j++)
        {
            std::cout<<(*this)(i,j);
            
            if (j < n-1)
            {
                std::cout<<",";
            }
        }
        std::cout<<"},";
        std::cout<<std::endl;
    }

}




void Vector::showVector(const char* str)
{
    std::cout<<"-------------------------------"<<std::endl;
    std::cout<<"VectorSize = "<<this->size<<std::endl;
    
    std::cout<<str<<"="<<std::endl;
    std::cout<<"{";
    for (unsigned int i = 0; i < this->size; i++)
    {
        std::cout<<(*this)(i);
            
        if (i < size-1)
        {
            std::cout<<",";
        }
    }
    std::cout<<"},";
    std::cout<<std::endl;

}




// This function is for teset
void Matrix::setMatrixDataTranspose(float datain[][3], const unsigned int rowSize, const unsigned int colSize)
{
    ///---------------------
    /// float datain[rowIndex][colIndex]
    //
    
    m = colSize;
    n = rowSize;
    
    allocate(colSize,rowSize);
    for (unsigned int i = 0 ;  i < m ; i++)
    {
        for (unsigned int j = 0 ; j < n ; j++)
        {
            data[i+j*m] = (double)datain[j][i];
        }
    }
}

#endif




void Matrix::setMatrixData(double* datain, const unsigned int rowSize, const unsigned int colSize)
{
    ///---------------------
    /// double datain[rowIndex+colIndex*rowSize]
    //
    m = rowSize;
    n = colSize;
    
    allocate(rowSize, colSize);
    memcpy(data,datain,m*n*sizeof(double));
}




void Matrix::setMatrixData(float** datain, const unsigned int rowSize, const unsigned int colSize)
{
    ///---------------------
    /// float datain[rowIndex][colIndex]
    //
    m = rowSize;
    n = colSize;
    
    allocate(rowSize, colSize);
    for (unsigned int i = 0 ;  i < m ; i++)
    {
        for (unsigned int j = 0 ; j < n ; j++)
        {
            data[i+j*m] = (double)datain[i][j];
        }
    }
}



void Matrix::setMatrixData(double** datain, const unsigned int rowSize, const unsigned int colSize)
{
    ///---------------------
    /// float datain[rowIndex][colIndex]
    //
    
    m = rowSize;
    n = colSize;
    
    allocate(rowSize, colSize);
    for (unsigned int i = 0 ;  i < m ; i++)
    {
        for (unsigned int j = 0 ; j < n ; j++)
        {
            data[i+j*m] = datain[i][j];
            //printf("data[i+j*m]= , ", data[i+j*m]);
        }
    }
}




void Matrix::setMatrixDataTranspose(float** datain, const unsigned int rowSize, const unsigned int colSize)
{
    ///---------------------
    /// float datain[rowIndex][colIndex]
    //
    
    m = colSize;
    n = rowSize;
    
    allocate(colSize,rowSize);
    for (unsigned int i = 0 ;  i < m ; i++)
    {
        for (unsigned int j = 0 ; j < n ; j++)
        {
            data[i+j*m] = (double)datain[j][i];
        }
    }
}




void Matrix::setMatrixDataTranspose(double** datain, const unsigned int rowSize, const unsigned int colSize)
{
    ///---------------------
    /// float datain[rowIndex][colIndex]
    //
    
    m = colSize;
    n = rowSize;
    
    allocate(colSize,rowSize);
    for (unsigned int i = 0 ;  i < m ; i++)
    {
        for (unsigned int j = 0 ; j < n ; j++)
        {
            data[i+j*m] = datain[j][i];
        }
    }
}





void Matrix::setColumnVectorDataAndRepeatN(float vecDataIn[], const unsigned int vecSize, const unsigned int repMatSize)
{
    m = vecSize;
    n = repMatSize;
    
    allocate(vecSize,repMatSize);
    
    
    for (unsigned int i = 0 ;  i < m ; i++)
    {
        for (unsigned int j = 0 ; j < n ; j++)
        {
            data[i+j*m] = (double)vecDataIn[i];
        }
    }
}



void Matrix::sliceMatRowData(double* datain, const unsigned int rowSizeIn,const unsigned int rowSizeOut, const unsigned int colSize)
{
    ///---------------------
    /// double datain[rowIndex+colIndex*rowSize]
    //  the output matrix taks first [rowSizeOut] rows  from the input matrix
    // - the input matrix has  [rowSizeIn by colSize]
    // - the output matrix has [rwoSizeOut by colSize]
    //--------------------
    
    m = rowSizeOut;
    n = colSize;
    
    allocate(rowSizeOut, colSize);
    for (unsigned int j=0;j < colSize;j++)
    {
        memcpy(&data[j*rowSizeOut],&datain[j*rowSizeIn],rowSizeOut*sizeof(double));
    }
}



void Matrix::sliceMatColData(double* datain, const unsigned int colSizeIn,const unsigned int colSizeOut, const unsigned int rowSize)
{
    ///---------------------
    /// double datain[rowIndex+colIndex*rowSize]
    //  the output matrix taks first [colSizeOut] cols  from the input matrix
    // - the input matrix has  [rowSize by colSizeIn]
    // - the output matrix has [rwoSize by colSizeOut]
    //--------------------
    
    m = rowSize;
    n = colSizeOut;
    
    allocate(rowSize, colSizeOut);
    for (unsigned int i=0;i < rowSize;i++)
    {
        for (unsigned int j =0 ; j < colSizeOut ; j++)
        {
            (*this)(i,j) = datain[i + j*rowSize];
        }
    }
}



// memory allocation
void Matrix::allocate(unsigned int m_, unsigned int n_)
{
    // if already allocated, memory is freed
    deallocate();
    
    // new sizes
    m = m_;
    n = n_;
    
    data = new double[m_*n_];
    memset(data,0.0,m_*n_*sizeof(double));
    //std::cout<<"[Matrix] Allocate [double: "<<m<<" x "<<n<<"] memory for matrix storage"<<std::endl;
    
} // allocate




// memory free
void Matrix::deallocate()
{
    
    if (data != nullptr)
    {
        delete [] data;
        data = nullptr;
    }
    //std::cout<<"[Matrix] deAllocate memory"<<std::endl;
    
}



# ifdef __DSPLIBTEST__
// L2-norm ||A-B||^2

double matrix_compare(const Matrix& A, const Matrix& B)
{
    // matrices must have same size
    if (A.m != B.m or  A.n != B.n)
        return std::numeric_limits<double>::max();
    
    double res=0;
    for(unsigned int i = 0; i < A.m; i++) {
        for (unsigned int j = 0; j < A.n; j++) {
            res += (A(i,j)-B(i,j)) * (A(i,j)-B(i,j));
        }
    }
    
    res /= A.m*A.n;
    return res;
}



double mseVectorCompare(Vector& inVec, Vector& refVec)
{
    double mse = 0.0;
    double l2normsqr = 0.0;
    
    if (inVec.size == refVec.size)
    {
        for (unsigned int i = 0 ; i < inVec.size ; i++)
        {
            mse         += pow(inVec(i) - refVec(i),2.0);
            l2normsqr   += pow(refVec(i),2.0);
        }
        
        LOGD("[Vector] mse = %1.2f dB",10*log10(mse/l2normsqr));
        return mse/l2normsqr;
    }
    else
    {
        LOGD("[Vector] Different Vector size");
        return -1.0;
    }
}

double mseMatrixCompare(Matrix& inMat, Matrix& refMat)
{
    double acc = 0.0;
    double mse = 0.0;
    double l2normsqr = 0.0;
    
    
    if ((inMat.m == refMat.m) && (inMat.n == refMat.n))
    {
        for (unsigned int i = 0 ; i < inMat.m ; i++)
        {
            for(unsigned int j =0 ; j < inMat.n ; j++)
            {
                acc         += pow(inMat(i,j) - refMat(i,j),2.0);
                l2normsqr   += pow(refMat(i,j),2.0);
            }
        }
        
        if (l2normsqr > 0.0)
        {
            mse = acc/l2normsqr;
            LOGD("[Matrix] Normalized mse = %1.2f dB",10*log10(mse));

        }
        else
        {
            mse = acc;
            LOGD("[Matrix] mse = %1.2f dB",10*log10(mse));

        }

        return mse;
    }
    else
    {
        LOGD("[Matrix] Different Matrix size");
        return -1.0;
    }
}
#endif


void householder(Matrix& mat, Matrix& R, Matrix& Q)
{
    unsigned int m = mat.m;
    unsigned int n = mat.n;
    
    // array of factor Q1, Q2, ... Qm
    //std::vector<Matrix> qv(m);
    
    Matrix* qv = new Matrix[m];
    
    // temp array
    Matrix z(mat);
    Matrix z1;
    
    for (int k = 0; k < n && k < m - 1; k++) {
        
        Vector e(m), x(m);
        double a;
        
        // compute minor
        z1.compute_minor(z, k);
        
        // extract k-th column into x
        z1.extract_column(x, k);
        
        a = x.norm();
        if (mat(k,k) > 0) a = -a;
        
        for (int i = 0; i < e.size; i++)
            e(i) = (i == k) ? 1 : 0;
        
        // e = x + a*e
        vmadd(x, e, a, e);
        
        // e = e / ||e||
        e.rescale_unit();
        
        // qv[k] = I - 2 *e*e^T
        compute_householder_factor(qv[k], e);
        
        // z = qv[k] * z1
        z.mult(qv[k], z1);
        
    }
    
    Q = qv[0];
    
    // after this loop, we will obtain Q (up to a transpose operation)
    for (int i = 1; i < n && i < m - 1; i++)
    {
        
        z1.mult(qv[i], Q);
        Q = z1;
        
    }
    
    R.mult(Q, mat);
    Q.transpose();
    delete [] qv;
}




/* sparse matrix multiplicatio:  A*B = C */

# ifdef __DSPLIBTEST__
bool mulSparseMatrix(SparseMatrix* Ain, SparseMatrix* Bin, SparseMatrix* Cout)
{
    // implementation of matrix multiplication of
    //  C = A*B
    //  where C is a [M by N] sparse matrix
    //        A is a [M by L] sparse matrix
    //        B is a [L by N] sparse matrix
    // This code is written by Jaewook Kang @ 2017 May
    
    //1 feasibility check!
    // feasibility check
    if (Ain->getColSize() !=  Bin->getRowSize() )
    {
        LOGE("[sparseMatrix] Invalid matrix muliplication!");
        return false;
    }
    
    unsigned int colSizeCout = Bin->getColSize();
    unsigned int rowSizeCout = Ain->getRowSize();


    SparseMatrixRepTable* tempTableForCout = nullptr;
    
    // temporal memory allocation
    tempTableForCout = new SparseMatrixRepTable[colSizeCout*rowSizeCout];
    memset(tempTableForCout,0.0,sizeof(SparseMatrixRepTable)*colSizeCout*rowSizeCout);
    
    
    // matrix multiplication
    unsigned int firstNonZeroIndexInARow = 0;
    unsigned int prevFirstNonZeroIndexInARow = 0;
    unsigned int nonZeroCntInArow = 0;
    unsigned int nonZeroCntForCout = 0;
    double* tempCoutRowVec =  new double[colSizeCout];

    
    // row loop
    while(true)
    {
        firstNonZeroIndexInARow = prevFirstNonZeroIndexInARow + nonZeroCntInArow;
        
        // m is Cout rowIndex
        int m  = Ain->getRowIndex(firstNonZeroIndexInARow);
        
        if (m == -1)
        {// loop termination condition
            break;
        }

        //----  loop for inner product of vectors
        nonZeroCntInArow = 0;
        memset(tempCoutRowVec,0.0,sizeof(double)*colSizeCout);

        while(Ain->getRowIndex(firstNonZeroIndexInARow + nonZeroCntInArow) == m)
        {
            unsigned int tempNonZeroIndexOfA = firstNonZeroIndexInARow + nonZeroCntInArow;
            unsigned int colIndexOfA = Ain->getColIndex(tempNonZeroIndexOfA);
            double tempElemOfA       = Ain->getElemValue(tempNonZeroIndexOfA);
            
            double tempElemOfB = 0.0;
            

            // column loop
            for(unsigned int n = 0 ; n < colSizeCout; n++) // n is Cout columIndex
            {
                tempElemOfB     = Bin->getElemValue(colIndexOfA,n);
                if (tempElemOfB!=0.0)
                {
                    tempCoutRowVec[n] = tempCoutRowVec[n] + tempElemOfA *tempElemOfB;
                }
  
            }
            nonZeroCntInArow++;
            
        }
        //------ inner product loop end
        for(unsigned int n = 0 ; n < colSizeCout; n++) // n is Cout columIndex
        {
            double elemCout = tempCoutRowVec[n];

            if (elemCout != 0.0)
            {
                tempTableForCout[nonZeroCntForCout].rowIndex = m;
                tempTableForCout[nonZeroCntForCout].colIndex = n;
                tempTableForCout[nonZeroCntForCout].values = elemCout;
                nonZeroCntForCout++;
            }
        }
        
        prevFirstNonZeroIndexInARow = firstNonZeroIndexInARow;
    }
    Cout->setMatrix(tempTableForCout,rowSizeCout,colSizeCout,nonZeroCntForCout);
    
    
    // memory deallocation
    delete [] tempTableForCout;
    delete [] tempCoutRowVec;
    tempTableForCout = nullptr;
    tempCoutRowVec   = nullptr;
    
    return true;
}
#endif




bool mulSparseMatrix(SparseMatrix* Ain, Matrix* Bin, Matrix* Cout)
{
    // implementation of matrix multiplication of
    //  C = A*B
    //  where C is a [M by N] dense matrix
    //        A is a [M by L] sparse matrix
    //        B is a [L by N] dense matrix
    // use an algorithm from https://www.cs.cmu.edu/~scandal/cacm/node9.html
    //
    // c_mn =  sum ( a_row_m[l] * b_col_n[l] :  a_row_m[l] is nonzero )
    // This code is written by Jaewook Kang @ 2017 Oct
    //

    // feasibility check
    if (Ain->getColSize() !=  Bin->getRowSize() )
    {

        LOGE("[sparseMatrix] Invalid matrix muliplication!");
        return false;
    }
    
    unsigned int colSizeCout = Bin->getColSize();
    unsigned int rowSizeCout = Ain->getRowSize();

    Cout->allocate(rowSizeCout, colSizeCout);

    // matrix multiplication
    unsigned int firstNonZeroIndexInARow = 0;
    unsigned int prevFirstNonZeroIndexInARow = 0;
    unsigned int nonZeroCntInArow = 0;
    
    // row loop
    while(true)
    {
        firstNonZeroIndexInARow = prevFirstNonZeroIndexInARow + nonZeroCntInArow;
        
        // m is Cout rowIndex
        int m  = Ain->getRowIndex(firstNonZeroIndexInARow);
        
        if (m == -1)
        {// loop termination condition
            break;
        }
        
        //----  loop for inner product of vectors
        nonZeroCntInArow = 0;
        while(Ain->getRowIndex(firstNonZeroIndexInARow + nonZeroCntInArow) == m)
        {
            unsigned int tempNonZeroIndexOfA = firstNonZeroIndexInARow + nonZeroCntInArow;
            unsigned int colIndexOfA = Ain->getColIndex(tempNonZeroIndexOfA);
        
            double tempElemOfA       = Ain->getElemValue(tempNonZeroIndexOfA);
            double tempElemOfB = 0.0;
            
            // column loop
            for(unsigned int n = 0 ; n < colSizeCout; n++) // n is Cout columIndex
            {
                tempElemOfB   = (*Bin)(colIndexOfA,n);
                if (tempElemOfB != 0.0)
                {
                    (*Cout)(m,n)     = (*Cout)(m,n) + tempElemOfA * tempElemOfB;
                }
            }
            nonZeroCntInArow++;

        }
        //------ inner product loop end
    
        prevFirstNonZeroIndexInARow = firstNonZeroIndexInARow;
    }
    return true;
}

bool isEqualMatSize(Matrix* Ain, Matrix* Bin)
{
    if ( (Ain->getColSize() == Bin->getColSize() ) && (Ain->getRowSize() == Bin->getRowSize()))
        return true;
    else
        return false;
    
}


bool elementWiseMatMul(Matrix* Ain, Matrix* Bin, Matrix* Cout)
{
    if ( isEqualMatSize(Ain, Bin) == false)
    {
        LOGE("[EverybodyMathLib] The two matrices have difference dimension.");
        return false;
    }
    else
    {
    
        unsigned int rowSize = Ain->getRowSize();
        unsigned int colSize = Ain->getColSize();
        
        for (unsigned int i = 0 ;  i < rowSize ; i++)
        {
            for (unsigned int j = 0 ; j < colSize ; j++)
            {
                (*Cout)(i,j) = (*Ain)(i,j) * (*Bin)(i,j);
            }
        }
        return true;
    }
}


bool elementWiseMatDiv(Matrix* Ain, Matrix* Bin, Matrix* Cout)
{
    if ( isEqualMatSize(Ain, Bin) == false)
    {
        LOGE("[EverybodyMathLib] The two matrices have difference dimension.");
        return false;
    }
    else{
    
        unsigned int rowSize = Ain->getRowSize();
        unsigned int colSize = Ain->getColSize();
        
        for (unsigned int i = 0 ;  i < rowSize ; i++)
        {
            for (unsigned int j = 0 ; j < colSize ; j++)
            {
                (*Cout)(i,j) = (*Ain)(i,j) / (*Bin)(i,j);

            }
        }
        return true;
    }
}

double innerProduct(Vector& Ain, Vector& Bin)
{
    unsigned int Asize = Ain.size;
    unsigned int Bsize = Bin.size;
    
    double out = 0.0;
    
    if (Asize != Bsize)
    {
        LOGD("[EverybodyMathLib] Two input vectors have different size.");
    }
    else
    {
        for (unsigned int i = 0 ; i <Asize ; i++)
        {
            out += Ain(i) * Bin(i);
        }
    }
    return out;

}

void diagMatrix(Vector& VecIn, Matrix& MatOut)
{
    unsigned int size = VecIn.size;
    MatOut.allocate(size, size);
    
    for (unsigned int i = 0 ;  i  < size ; i++ )
    {
        MatOut(i,i) = VecIn(i);
    }
    
}

