//
//  test_runner.cpp
//  MathLib
//
//  Created by jwkangmacpro on 2018. 5. 23..
//  Copyright © 2018년 jwkangmacpro. All rights reserved.
//
#include "gtest/gtest.h"

#include "header/EverybodyDSPLib.h"

#define TEST_SPARSE_MTX_MUL
#define TEST_MTX_MUL
#define TEST_IMPORT_MTX
#define TEST_QRDECMP
#define TEST_SVD



#ifdef TEST_SPARSE_MTX_MUL

TEST(uTEST_sparsematrix, test1)
{
    
    const int rowSize0 = 4;
    const int colSize0 = 5;
    const char* intputfilename = "./data_for_sparseMatrix/matrix_input.txt";
    
    SparseMatrix mat0(rowSize0,colSize0);
    mat0.readMatrix(intputfilename);
    mat0.showMatrix("mat0");
    
    LOGD("[SparseMatrix] Transpose\n");
    mat0.transpose();
    mat0.showMatrix("mat0");
    
    SparseMatrix mat1;
    mat1.transpose(&mat0);
    mat1.showMatrix("mat1");
}

TEST(uTEST_sparsematrix, test2)
{
    const int rowSize0 = 20;
    const int colSize0 = 10;
    const char* intputfilename = "./data_for_sparseMatrix/matrix_input2.txt";
    
    SparseMatrix mat0(rowSize0,colSize0);
    mat0.readMatrix(intputfilename);
    mat0.showMatrix("mat0");
    
    LOGD("[SparseMatrix] Transpose\n");
    //    mat0.transpose();
    //    mat0.showMatrix("mat0");
}

TEST(uTEST_sparsematrix, test3)
{
    const int rowSize0 = 20;
    const int colSize0 = 20;
    const char* intputfilename = "./data_for_sparseMatrix/matrix_input3.txt";
    
    SparseMatrix mat0(rowSize0,colSize0);
    mat0.readMatrix(intputfilename);
    mat0.showMatrix("mat0");
    
    //    LOGD("[SparseMatrix] Transpose\n");
    //    mat0.transpose();
    //    mat0.showMatrix("mat0");
}



TEST(uTEST_sparsematrix, sparse_sparse_mul)
{
    const int rowSize0 = 4;
    const int colSize0 = 5;
    const char* intputfilename = "./data_for_sparseMatrix/matrix_input.txt";
    
    const int rowSize1 = 5;
    const int colSize1 = 4;
    const char* intputfilename1 = "./data_for_sparseMatrix/matrix_input4.txt";
    LOGD("=============================================================\n");
    LOGD("# [Main] Reading the first Matrix, denoted by A, from a file\n");
    
    SparseMatrix mat0(rowSize0,colSize0);
    mat0.readMatrix(intputfilename);
    mat0.showMatrix("mat0");
    
    LOGD("=============================================================\n");
    LOGD("# [Main] Reading the second Matrix, denoted by B, from a file\n");
    
    SparseMatrix mat1(rowSize1,colSize1);
    mat1.readMatrix(intputfilename1);
    //mat1.showMatrix("mat1");
    
    SparseMatrix matOut(rowSize0,colSize1);;
    LOGD("=============================================================\n");
    LOGD("# [Main] Sparse X Sparse Matrix mutiplications, denoted by C=A*B\n");
    //    matOut =  mat0 * mat1
    mulSparseMatrix(&mat0, &mat1, &matOut);
    matOut.showMatrix("matOut");
    /* correct ans
     matOut =
     20     6     0     0
     0    10    28     4
     0     0     0     0
     2    12     0    12
     */
    
}

TEST(uTEST_sparsematrix, sparse_dense_mul)
{
    const int rowSize0 = 4;
    const int colSize0 = 5;
    const char* intputfilename = "./data_for_sparseMatrix/matrix_input.txt";
    
    const int rowSize1 = 5;
    const int colSize1 = 4;
    const char* intputfilename1 = "./data_for_sparseMatrix/matrix_input4.txt";
    LOGD("=============================================================\n");
    LOGD("# [Main] Reading the first Matrix, denoted by A, from a file\n");
    
    SparseMatrix mat0(rowSize0,colSize0);
    mat0.readMatrix(intputfilename);
    //    mat0.showMatrix("mat0");
    
    LOGD("=============================================================\n");
    LOGD("# [Main] Reading the second Matrix, denoted by B, from a file\n");
    
    Matrix mat1(rowSize1,colSize1);
    mat1.readMatrix(intputfilename1);
    //mat1.showMatrix("mat1");
    
    Matrix matOut(rowSize0,colSize1);
    LOGD("=============================================================\n");
    LOGD("# [Main] Sparse X dense Matrix mutiplications, denoted by C=A*B\n");
    //    matOut =  mat0 * mat1
    
    mulSparseMatrix(&mat0, &mat1, &matOut);
    matOut.showMatrix("matOut");
    /* correct ans
     matOut =
     20     6     0     0
     0    10    28     4
     0     0     0     0
     2    12     0    12
     */
}
#endif //TEST_SPARSE_MTX_MUL

#ifdef TEST_MTX_MUL

TEST(uTEST_matrix_elementwiseMulDiv,test1)
{
    double Aindata[][3] = {
        { 12, -51,   4},
        {  6, 167, -68},
        { -4,  24, -41},
        { -4,  24, -41},
        
    };
    
    double Bindata[][3] = {
        { 1, 0,   4},
        {  6, 2, 0},
        { -4,  5, 3},
    };
    
    Matrix Ain(Aindata);
    Matrix Aint;
    Ain.showMatrix("Ain");
    //    Aint.transpose(Ain);
    //    Aint.showMatrix("Aint");
    Ain.transpose();
    Ain.showMatrix("Ain");
    std::cout<<"-----------------------------"<<std::endl;
    Matrix Bin(Bindata);
    Bin.showMatrix("Bin");
    std::cout<<"-----------------------------"<<std::endl;
    
    Matrix CMulOut,CDivOut;
    
    elementWiseMatDiv(&Ain, &Bin, &CDivOut);
    elementWiseMatMul(&Ain, &Bin, &CMulOut);
    
    CMulOut.showMatrix("CMulOut");
    std::cout<<"-----------------------------"<<std::endl;
    CDivOut.showMatrix("CDivOut");
    
    
    
}


TEST(uTEST_matrix_setMatrixDataTranspose, test1)
{
    float Xmix[][3] = {
        { 12, -51,   4},
        {  6, 167, -68},
        { -4,  24, -41},
        { -1,   1,   0},
        {  2,   0,   3},
    };
    
    Matrix X;
    X.setMatrixDataTranspose(Xmix,5,3);
    X.showMatrix("X");
    
    
}


TEST(uTEST_matrix_setColumnVectorDataAndRepeatN,test1)
{
    float vec[5] = {1,2,3,4,5};
    
    Matrix X;
    X.setColumnVectorDataAndRepeatN(vec, 5, 3);
    X.showMatrix("X");
}
#endif //TEST_MTX_MUL


#ifdef TEST_QRDECMP


TEST(everybodyMathLibTEST, qrdmp)
{
    
    double in[][3] = {
        { 12, -51,   4},
        {  6, 167, -68},
        { -4,  24, -41},
        { -1,   1,   0},
        {  2,   0,   3},
    };
    
    Matrix A(in);
    Matrix Q, R;
    
    A.showMatrix("A");
    
    // compute QR decompostion
    householder(A, R, Q);
    
    Q.showMatrix("Q");
    R.showMatrix("R");
    
    // compare Q*R to the original matrix A
    Matrix A_check;
    A_check.mult(Q, R);
    
    // compute L2 norm ||A-A_check||^2
    //    double l2 = matrix_compare(A,A_check);
    
    // display Q*R
    //    A_check.showMatrix(l2 < 1e-12 ? "A == Q * R ? yes" : "A == Q * R ? no");
    
}
#endif //TEST_QRDECMP



#ifdef TEST_SVD

TEST(jacobi_svd,test1)
{
    unsigned int rowsize = 11;
    unsigned int colsize = 5;
    bool matrixExportOptOn = false;
    bool qropt = false;
    
    const char* inputfilename1 = "./data_for_utest_jacobisvd/inputmatrixA_test1.txt";
    const char* reffilename1   = "./data_for_utest_jacobisvd/matrixU_test1.txt";
    const char* reffilename2   = "./data_for_utest_jacobisvd/matrixV_test1.txt";
    const char* reffilename3   = "./data_for_utest_jacobisvd/singularvalue_test1.txt";
    
    
    Vector singularValuemat(colsize);
    singularValuemat.readVector(reffilename3);
    
    Matrix Ain(rowsize,colsize);
    Ain.readMatrix(inputfilename1);
    
    Matrix Umat(rowsize,colsize);
    Umat.readMatrix(reffilename1);
    
    Matrix Vmat(colsize,colsize);
    Vmat.readMatrix(reffilename2);
    
    JacobiSVD svdworker(Ain,matrixExportOptOn);
    
    svdworker.showMatrixA();
    svdworker.compute(qropt);
    //svdworker.showResult();
    
    LOGD("[TESTRUNNER] Singular values MSE dB");
    mseVectorCompare(svdworker.getSingularValues(),singularValuemat);
    
    if (matrixExportOptOn)
    {
        LOGD("[TESTRUNNER] Matrix U MSE dB");
        mseMatrixCompare(svdworker.getMatrixU(), Umat);
        LOGD("[TESTRUNNER] Matrix V MSE dB");
        mseMatrixCompare(svdworker.getMatrixV(), Vmat);
    }
}


TEST(jacobi_svd,test2)
{
    unsigned int rowsize = 25;
    unsigned int colsize = 5;
    bool matrixExportOptOn = false;
    bool qropt = false;
    
    const char* inputfilename1 = "./data_for_utest_jacobisvd/inputmatrixA_test2.txt";
    const char* reffilename1   = "./data_for_utest_jacobisvd/matrixU_test2.txt";
    const char* reffilename2   = "./data_for_utest_jacobisvd/matrixV_test2.txt";
    const char* reffilename3   = "./data_for_utest_jacobisvd/singularvalue_test2.txt";
    
    
    Vector singularValuemat(colsize);
    singularValuemat.readVector(reffilename3);
    
    Matrix Ain(rowsize,colsize);
    Ain.readMatrix(inputfilename1);
    
    Matrix Umat(rowsize,colsize);
    Umat.readMatrix(reffilename1);
    
    Matrix Vmat(colsize,colsize);
    Vmat.readMatrix(reffilename2);
    
    JacobiSVD svdworker(Ain,matrixExportOptOn);
    
    svdworker.showMatrixA();
    svdworker.compute(qropt);
    //svdworker.showResult();
    
    LOGD("[TESTRUNNER] Singular values MSE dB");
    mseVectorCompare(svdworker.getSingularValues(),singularValuemat);
    
    if (matrixExportOptOn)
    {
        LOGD("[TESTRUNNER] Matrix U MSE dB");
        mseMatrixCompare(svdworker.getMatrixU(), Umat);
        LOGD("[TESTRUNNER] Matrix V MSE dB");
        mseMatrixCompare(svdworker.getMatrixV(), Vmat);
    }
}


TEST(jacobi_svd,test5)
{
    unsigned int rowsize = 25;
    unsigned int colsize = 5;
    bool matrixExportOptOn = false;
    bool qropt = false;
    
    const char* inputfilename1 = "./data_for_utest_jacobisvd/inputmatrixA_test5.txt";
    const char* reffilename1   = "./data_for_utest_jacobisvd/matrixU_test5.txt";
    const char* reffilename2   = "./data_for_utest_jacobisvd/matrixV_test5.txt";
    const char* reffilename3   = "./data_for_utest_jacobisvd/singularvalue_test5.txt";
    
    
    Vector singularValuemat(colsize);
    singularValuemat.readVector(reffilename3);
    
    Matrix Ain(rowsize,colsize);
    Ain.readMatrix(inputfilename1);
    
    Matrix Umat(rowsize,colsize);
    Umat.readMatrix(reffilename1);
    
    Matrix Vmat(colsize,colsize);
    Vmat.readMatrix(reffilename2);
    
    JacobiSVD svdworker(Ain,matrixExportOptOn);
    
    svdworker.showMatrixA();
    svdworker.compute(qropt);
    svdworker.showResult();
    
    LOGD("[TESTRUNNER] Singular values MSE dB");
    mseVectorCompare(svdworker.getSingularValues(),singularValuemat);
    
    if (matrixExportOptOn)
    {
        LOGD("[TESTRUNNER] Matrix U MSE dB");
        mseMatrixCompare(svdworker.getMatrixU(), Umat);
        LOGD("[TESTRUNNER] Matrix V MSE dB");
        mseMatrixCompare(svdworker.getMatrixV(), Vmat);
    }
}

#endif // TEST_SVD

int main(int argc, char* argv[])
{
    printf("Successful Running in Google Test\n");
    
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
    
    return 0;
}


