//
//  EverybodyMathLib.h
//
//  Created by jwkangmacpro on 2017. 4. 18..
//
//  The first reference :  https://rosettacode.org/wiki/QR_decomposition

// update history
//  The first revision from the reference is done by jwkang 2017 Apr.
//  - addition of matrix algebra @ 2018 Fev
//
//  final update @ 2018 feb by jaewookkang



#ifndef EverybodyMathLib_h
#define EverybodyMathLib_h

#include "log.h"
#include <cmath>
#include <cstring> // for memset
#include <cstdio>
//#include <memory>
#include "sparseMatrix.h"

# ifdef __DSPLIBTEST__
#include <fstream>// not used for android
#include <iostream>
#endif


const long double PI=3.1415926535897932;


template<class InputIterator, class OutputIterator>
OutputIterator Copy_v2(InputIterator first, InputIterator last, OutputIterator result)
{
    while (first!=last) {
        *result = *first;
        ++result; ++first;
    }
    return result;
}


// column vector
class Vector {
    
public:
    // default constructor (don't allocate)
    Vector() : size(0), data(nullptr) {}
    
    // constructor with memory allocation, initialized to zero
    Vector(int size_) : Vector() {
        size = size_;
        allocate(size_);
    }
    
    // destructor
    ~Vector() {
        deallocate();
    }
    
    // access data operators
    double& operator() (int i) {
        return data[i]; }
    double  operator() (int i) const {
        return data[i]; }
    
    // operator assignment
    Vector& operator=(const Vector& source) {
        
        // self-ass&ignment check
        if (this != &source) {
            if ( size != (source.size) ) {   // storage cannot be reused
                allocate(source.size);         // re-allocate storage
            }
            // storage can be used, copy data
//            std::copy(source.data, source.data + source.size, data);
            Copy_v2(source.data, source.data + source.size, data);

        }
        return *this;
    }
    
    // memory allocation
    void allocate(int size_) {
        
        deallocate();
        
        // new sizes
        size = size_;
        
        data = new double[size_];
        memset(data,0,size_*sizeof(double));
        
    } // allocate
    
    // memory free
    void deallocate() {
        
        if (data != nullptr)
        {
            delete [] data;
            data = nullptr;
        }

    }
    
    //   ||x||
    double norm() {
        double sum = 0;
        for (int i = 0; i < size; i++) sum += (*this)(i) * (*this)(i);
        return sqrt(sum);
    }
    
    // divide data by factor
    void rescale(double factor) {
        for (int i = 0; i < size; i++) (*this)(i) /= factor;
    }
    
    void rescale_unit() {
        double factor = norm();
        rescale(factor);
    }
    
    void setData(double* datain)
    {
        for (unsigned int i = 0; i < size ; i++)
        {
            (*this)(i) =  datain[i];
        }
        
    }
    
    void powElem(const double exponent);    // pow(v)
    void absElem();                         // abs(v)
    void absPowElem(const double exponent); // pow(abs(v))
    void absSqrtElem();                     // sqrt(abs(v))
    
    double sumAllElem();                            // sum(v)
    double powSumAllElem(const double exponent);    // sum(pow(v,exponent))
    double absPowSumAllElem(const double exponent); // sum(pow(abs(v),exponent))
    
    
# ifdef __DSPLIBTEST__
    void showVector(const char* str);
    void readVector(const char* inputfilename);
#endif
    
    int size;
    
private:
    double* data;
    
}; // class Vector



class Matrix {
    
public:
    // default constructor (don't allocate)
    Matrix() : m(0), n(0), data(nullptr) {}
    
    // constructor with memory allocation, initialized to zero
    Matrix(int m_, int n_) : Matrix()
    {
        m = m_; // rowSize
        n = n_; // colSize
        allocate(m_,n_);
    }
    
    // copy constructor
    Matrix(const Matrix& mat) : Matrix(mat.m,mat.n)
    {
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                (*this)(i,j) = mat(i,j);
    }
    
    // constructor from double array
    template<int rows, int cols>
    Matrix(double (&a)[rows][cols]) : Matrix(rows,cols) {
        
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                (*this)(i,j) = a[i][j];
    }
    
    // constructor from float array
    template<int rows, int cols>
    Matrix(float (&a)[rows][cols]) : Matrix(rows,cols) {
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                (*this)(i,j) = (double)a[i][j];
    }
    
    // destructor
    ~Matrix() {
        deallocate();
    }
    
    
    // access data operators
    double& operator() (unsigned int i, unsigned int j) {
        return data[i+m*j]; }
    double  operator() (unsigned int i, unsigned int j) const {
        return data[i+m*j]; }
    
    // operator assignment
    Matrix& operator=(const Matrix& source) {
        
        // self-assignment check
        if (this != &source) {
            if ( (m*n) != (source.m * source.n) ) { // storage cannot be reused
                allocate(source.m,source.n);          // re-allocate storage
            }
            // storage can be used, copy data
//            std::copy(source.data, source.data + source.m*source.n, data);
            Copy_v2(source.data, source.data + source.m*source.n, data);
        }
        return *this;
    }
    
    // compute minor
    void compute_minor(const Matrix& mat, unsigned int d) {
        
        allocate(mat.m, mat.n);
        
        for (unsigned int i = 0; i < d; i++)
            (*this)(i,i) = 1.0;
        for (int i = d; i < mat.m; i++)
            for (unsigned int j = d; j < mat.n; j++)
                (*this)(i,j) = mat(i,j);
        
    }
    
    // Matrix multiplication
    // c = a * b
    // c will be re-allocated here
    void mult(const Matrix& a, const Matrix& b) {
        
        if (a.n != b.m) {
            LOGD("[EverybodyMathLib] Matrix multiplication not possible, sizes don't match !");
            return;
        }
        
        // reallocate ourself if necessary i.e. current Matrix has not valid sizes
        if (a.m != m or b.n != n)
            allocate(a.m, b.n);
        
        memset(data,0,m*n*sizeof(double));
        
        for (unsigned int i = 0; i < a.m; i++)
            for (unsigned int j = 0; j < b.n; j++)
                for (unsigned int k = 0; k < a.n; k++)
                    (*this)(i,j) += a(i,k) * b(k,j);
        
    }
    
    // Matrix sum
    // c = a + b
    void sum(const Matrix& a, const Matrix& b){
        if ( (a.n != b.n) || (a.m != b.m))
        {
            LOGD("[EverybodyMathLib] Matrix sum not possible, sizes don't match !");
            return;
        }
        
        for (unsigned int i = 0; i < m ; i++)
        {
            for(unsigned int j =0 ; j < n ; j++)
            {
                (*this)(i,j) = a(i,j) + b(i,j);
            }
        }
    }
    
    // Matrix weighted sum
    // c = (1-weighted) * a + weight * b
    void weightedSum(const Matrix& a, const Matrix& b, float weight){
        if ( (a.n != b.n) || (a.m != b.m))
        {
            LOGD("[EverybodyMathLib] Matrix sum not possible, sizes don't match !");
            return;
        }
        
        for (unsigned int i = 0; i < m ; i++)
        {
            for(unsigned  int j =0 ; j < n ; j++)
            {
                (*this)(i,j) = (1-weight)*a(i,j) + weight*b(i,j);
            }
        }
    }
    
    // matrix transpose
    // matrix is re-allocated here
    Matrix& transpose() {
        Matrix temp;
        temp.setMatrixData(data, m, n);
        allocate(n, m);// reset
        for (unsigned int i = 0; i < m; i++) {
            for (unsigned int j = 0; j < n; j++) {
                (*this)(i,j) = temp(j,i);
            }
        }
        
        return *this;
    }
    
    Matrix& transpose(const Matrix& a) {
        allocate(a.n, a.m);
        for (unsigned int i = 0; i < m; i++) {
            for (unsigned int j = 0; j < n ; j++) {
                (*this)(i,j) = a(j,i);
            }
        }
        return *this;
    }
    
    void eye(const int size)
    {
        m = size;
        n = size;
        allocate(m,n);
        for (int i = 0 ; i < m ; i++)
        {
            (*this)(i,i) = 1.0;
        }
        
    }
    
    void triu(Matrix& Ain)
    {
        allocate(Ain.getRowSize(), Ain.getColSize());
        for ( int i = 0 ; i < m - 1 ; i++)
        {
            for ( int j = i + 1 ; j < n ; j++)
            {
                (*this)(i,j) = Ain(i,j);
            }
        }
        
    }
    
    void diag(Vector& diagVec)
    {
        int vecsize =0;
        if ( m <= n)
        {
            vecsize = m;
        }
        else
        {
            vecsize = n;
        }
        diagVec.allocate(vecsize);
        
        for ( int i = 0; i < diagVec.size ; i++)
        {
            diagVec(i) = (*this)(i,i);
        }
        
    }

    
# ifdef __DSPLIBTEST__
    // These functions are for test--------
    void setMatrixDataTranspose(float datain[][3], const unsigned int rowSize, const unsigned int colSize);
#endif
    
    void setMatrixData(double* datain, const unsigned int rowSize, const unsigned int colSize);// not used
    void setMatrixData(double** datain, const unsigned int rowSize, const unsigned int colSize);
    void setMatrixData(float** datain, const unsigned int rowSize, const unsigned int colSize);

    void setMatrixDataTranspose(float** datain, const unsigned int rowSize, const unsigned int colSize);
    void setMatrixDataTranspose(double** datain, const unsigned int rowSize, const unsigned int colSize);
    
    void setColumnVectorDataAndRepeatN(float vecDataIn[], const unsigned int vecSize, const unsigned int repMatSize);
    
    // slicing matrix in row
    void sliceMatRowData(double* datain, const unsigned int rowSizeIn, const unsigned int rowSizeOut, const unsigned int colSize);
    void sliceMatColData(double* datain, const unsigned int colSizeIn, const unsigned int colSizeOut, const unsigned int rowSize);
    double* getMatrixDataInDoubleArray() const {return data;}

    // take c-th column of m, put in v
    void extract_column(Vector& v_extracted, unsigned int c);
    void insert_column(Vector& v_inserted, unsigned int c);
    
    // memory allocation
    void allocate(unsigned int m_, unsigned int n_);

    // memory free
    void deallocate();

    
    unsigned int m; // the number of rows
    unsigned int n; // the number of columns
    
    unsigned int getColSize() const { return n;}
    unsigned int getRowSize() const { return m;}
    
    void powElem(double exponent);  // pow(mat,exponent)
    void absElem();                 // abs(mat)
    void absSqrtElem();             // sqrt(abs(mat)
    
    double sumAllElem();            // sum(sum(mat))
    double powSumAllElem(const double exponent);    // sum(sum(pow(mat,exponent)))
    double absSumAllElem();                         // sum(sum(abs(mat)))
    
    void columnMulScaling(const unsigned int colindex, const double scalingConst)
    {
        if ( colindex <= this->getColSize())
        {
            for (int i = 0 ; i < m ; i++)
            {
                (*this)(i,colindex) = (*this)(i,colindex) * scalingConst;
            }
        }
        else
        {
            LOGD("[Matrix] ERROR. The colindex > colsize of Matrix!!");
        }
    }
    
    void columnDivScaling(const unsigned int colindex, const double scalingConst)
    {
        if ( colindex <= this->getColSize())
        {
            if ( fabs(scalingConst) > 0)
            {
                for (int i = 0 ; i < m ; i++)
                {
                    (*this)(i,colindex) = (*this)(i,colindex) / scalingConst;
                }
            }
            else
            {
                LOGD("[Matrix] Divide by zero problem occur!");
            }
        }
        else
        {
            LOGD("[Matrix] ERROR. The colindex > colsize of Matrix!!");
        }
    }
    

# ifdef __DSPLIBTEST__
    // These functions are for test--------
    void readMatrix(const char* inputfilename);
    void showMatrix(const char* str);

#endif
    
private:
    double* data;
    
}; // struct Matrix

// some operation definition
// c = a + b * s
void vmadd(const Vector& a, const Vector& b, double s, Vector& c);

// vecout = c1*a + c2*b
void vmadd(const Vector& a, const Vector& b, const double c1, const double c2, Vector& vecout);

// mat = I - 2*v*v^T
// !!! m is allocated here !!!
void compute_householder_factor(Matrix& mat, const Vector& v);

void householder(Matrix& mat, Matrix& R,Matrix& Q);


/* sparse matrix multiplicatio:  A*B = C */
bool mulSparseMatrix(SparseMatrix* Ain, Matrix* Bin, Matrix* Cout);

/* elementwise matrix operations */
bool elementWiseMatMul(Matrix* Ain, Matrix* Bin, Matrix* Cout);
bool elementWiseMatDiv(Matrix* Ain, Matrix* Bin, Matrix* Cout);
bool isEqualMatSize(Matrix* Ain, Matrix* Bin);

void diagMatrix(Vector& VecIn, Matrix& MatOut);


/* vector inner product */
double innerProduct(Vector& Ain, Vector& Bin);

# ifdef __DSPLIBTEST__
// These functions are for test--------
bool mulSparseMatrix(SparseMatrix* Ain, SparseMatrix* Bin, SparseMatrix* Cout);// not used
// L2-norm ||A-B||^2
double matrix_compare(const Matrix& A, const Matrix& B);// not used
double mseVectorCompare(Vector& inVec, Vector& refVec);
double mseMatrixCompare(Matrix& inMat, Matrix& refMat);

#endif



#endif /* everybodyMathLib_h */
