//
//  JacobiSvd.cpp
//
//  Created by jwkangmacpro on 2018. 2. 2..
//  Copyright © 2018년 jwkangmacpro. All rights reserved.
//

#include "JacobiSvd.h"


JacobiSVD::JacobiSVD()
{
    mColSize_A = 0;
    mRowSize_A = 0;
    mNumOfSingularValues = 0;
    mNuclearNorm         = 0.0;
    mSVDCompute_ms       = 0.0;
    //mRankA                 = 0;
    mIsSingularVectorMatrixExport   = false;
    mIsSvdRun                       = false;
}




JacobiSVD::JacobiSVD(const Matrix& Ain, const bool matrixExportOptOn)
{
    reset(Ain,matrixExportOptOn);
    printf("[JacobiSVD] call JacobiSVD.");

}




JacobiSVD::JacobiSVD(const int colsize, const int rowsize, double** Adata, const bool matrixExportOptOn)
{
    reset(colsize, rowsize, Adata,matrixExportOptOn);
}



JacobiSVD::~JacobiSVD()
{
    if (mIsSingularVectorMatrixExport)
    {
        mS.deallocate();
        mU.deallocate();
        mV.deallocate();
    }
    mA.deallocate();
    mSingularValues.deallocate();
    //LOGD("[JacobiSVD] release JacobiSVD.")

}




void JacobiSVD::reset(const Matrix& Ain, const bool matrixExportOptOn)
{
    mNumOfSingularValues = 0;
    mNuclearNorm         = 0.0;
    mSVDCompute_ms       = 0.0;

    //mRankA                 = 0;
    mIsSvdRun = false;
    
    setMatrixA(Ain);
    mIsSingularVectorMatrixExport = matrixExportOptOn;

    printf("[JacobiSVD] Class is reset to perform svd.");
}




void JacobiSVD::reset(const int colsize, const int rowsize, double** Adata, const bool matrixExportOptOn)
{
    mNumOfSingularValues = 0;
    mNuclearNorm         = 0.0;
    mSVDCompute_ms       = 0.0;

    //mRankA                 = 0;

    mIsSvdRun            = false;
    
    setMatrixA(colsize, rowsize, Adata);
    mIsSingularVectorMatrixExport = matrixExportOptOn;

    //LOGD("[JacobiSVD] Class is reset to perform svd.");

}




void JacobiSVD::setMatrixA(const Matrix& Ain)
{
    mA = Ain;
    if (Ain.getRowSize() < Ain.getColSize())
    {
        mA.transpose(mA);
    }

    mColSize_A  = Ain.getColSize();
    mRowSize_A  = Ain.getRowSize();
    //LOGD("[JacobiSVD] mRowSize_A of A= %d",mRowSize_A);
    //LOGD("[JacobiSVD] mColSize_A of A= %d",mColSize_A);
}




void JacobiSVD::setMatrixA(const int colsize, const int rowsize, double** data)
{
    //LOGD("[JacobiSVD] allocate Matrix A");
    if (rowsize < colsize)
    {
        mA.setMatrixDataTranspose(data,rowsize,colsize);
    }
    else
    {
        mA.setMatrixData(data, rowsize, colsize);
    }
    
    mColSize_A = mA.getColSize();
    mRowSize_A = mA.getRowSize();
    
    //LOGD("[JacobiSVD] mRowSize_A of A= %d",mRowSize_A);
    //LOGD("[JacobiSVD] mColSize_A of A= %d",mColSize_A);

}




bool JacobiSVD::compute(const bool qropt)
{
    /* -----------------------------------------------------------------------------
     * The singular value decomposition computed by this function is the same as
     * that computed by the function svd.m in MATLAB.
     * Give an input matrix A We try to find S,V,U matrice satisfying
     *
     *   A = U * S * V^T
     *
     * where 
     * - A is a [ M by N ] input matrix (M>=N)
     * - U is a [ M by N ] unitary matrix (M>=N)
     * - S is a [ N by N ] singular value matrix
     * - V is a [ N by N ] unitary matrix
     *----------------------------------------------------
     * This code is based on a Jacobi rotation algorithm, is more accurate than the Bidiagonalization SVD.
     * However, it is also slower.

     * In our svd algorithm, we have applied
     * 1) householder qr decomposition to reduce input size of the Jacobi algorithm
     * 2) zero matrix exception routine to avoid uneccessary computation
     *
     * Note that our JacobiSVD class only supports svd calculation for tall matrix input A.
     * Namely, the size of input matrix A must satisfy rowsize >= colsize
     *
     * written by Jaewook Kang @ 2017 2/9 
     * ------------------------------------------------------------------------
     */
    
    time_t startTime, endTime;
    startTime = clock();
 
    const float svdIterationTol = 1E-4;
    const unsigned int maxIter = 10;
    unsigned int M = mRowSize_A;
    unsigned int N = mColSize_A;
    unsigned int numIter = 0;
    
    double OnDiagCost   = 0.0;
    // We test on each sweep to make sure the Off diagonal
    // sum is reducing. If it does not reduce we stop. We
    // use 9999 as the initial value.
    double OffDiagCost      = 0.0;
    double prevOffDiagCost  = 9999.0;
    
    double b_rr =0.0, b_cc = 0.0, b_rc = 0.0;
    double abs_b_rc = 0.0;
    
    double tau = 0.0;
    double theta = 0.0, alpha=0.0, beta = 0.0;
    double signTau = 0.0;
    
    unsigned int Rcnt = 0;

    /* for QR decomp */
    Matrix tempQ; // M by M
    Matrix Q; // M by K
    Matrix R; // M by N
    

    /* for matrix update in jacobi rotation iteration */
    Matrix currU;   // M by N
    Matrix nextU;   // M by N
    Matrix Ut;      // N by M
    
    Matrix currV;   // N by N
    Matrix G;       // N by N Jacobi rotation matrix
    Matrix B;       // B= U^T * U (N by N)
    Matrix triuB;   // N by N
    
    Vector colVectorOfU_r; // N by 1
    Vector colVectorOfU_c; // N by 1
    Vector tempVector;
    
    if(!mIsSvdRun)
    {
        mIsSvdRun = true;
        if (qropt)
        {
            // with qropt == true, M == N
            M=N;
            //**** step 1. QR decomposition
            //LOGD("[JacobiSVD] Use QR decomp for SVD");
            
            // mA = Q*R
            householder(mA,R,tempQ);
            //tempQ.showMatrix("Q");
            //R.showMatrix("R");
            
            // U = R(:,1:M)
            nextU.sliceMatRowData(R.getMatrixDataInDoubleArray(), mRowSize_A, M,N);// where (M ==N)
        }
        else
        {
            nextU = mA;
        }
        //nextU.showMatrix("U"); // M by N

        // On = sum(sum(U.^2)) ./N;
        OnDiagCost = nextU.powSumAllElem(2.0)/ N ;
        mV.eye(N);
        //**** step 2. Jacobi rotation algorithm
        for (numIter = 0 ; numIter < maxIter; numIter++)
        {
            //LOGD("----------------------------------");
            //LOGD("[JacobiSVD] numIter= %d ",numIter);
            //LOGD("[JacobiSVD] OnDiagCost = %1.4f",OnDiagCost);
            //LOGD("[JacobiSVD] OffDiagCost = %1.4f",OffDiagCost);
            //LOGD("[JacobiSVD] prevOffDiagCost = %1.4f",prevOffDiagCost);
            
            if (OnDiagCost <= 0.0)
            {
                printf("[JacobiSVD] A is a zero matrix");
                break;
            }
            
            // We count the rotations so that we know if we have not done any
            // during a whole sweep.
            Rcnt = 0 ;

            for( unsigned int r = 0 ; r < N - 1 ; r++)
            {
                for (unsigned int c = r + 1 ; c < N ; c++)
                {
                    //LOGD("----------------------------------");
                    //LOGD(" (r,c)=(%d,%d)",r,c);
                    
                    currV = mV;
                    currU = nextU;
                    //currU.showMatrix("currU");
                    //currV.showMatrix("currV");
                    
                    // Calculate the three elements of the implicit matrix B that are needed to calculate a Jacobi rotation. Since B is Hermitian, the fourth element (b_cr) is not needed.
                    
                    /*
                     // MATLAB codes
                     b_rr = sum(abs(U(:,r)).^2); % Real value.
                     b_cc = sum(abs(U(:,c)).^2); % Real value.
                     b_rc = U(:,r)' * U(:,c);    % Same type as U.
                     */
                    
                    // extract column of matrix U
                    currU.extract_column(colVectorOfU_r, r);
                    currU.extract_column(colVectorOfU_c, c);
                    
                    //colVectorOfU_r.showVector("U(:,r)");
                    //colVectorOfU_c.showVector("U(:,c)");
                    
                    //inner product
                    b_rc = innerProduct(colVectorOfU_c, colVectorOfU_r);
                    
                    /*
                     // deprecated at 180212 by jwkang
                    colVectorOfU_r.absPowElem(2.0);
                    colVectorOfU_c.absPowElem(2.0);
                    b_rr = colVectorOfU_r.sumAllElem();
                    b_cc = colVectorOfU_c.sumAllElem();
                    */
                    
                    b_rr = colVectorOfU_r.absPowSumAllElem(2.0);
                    b_cc = colVectorOfU_c.absPowSumAllElem(2.0);
                    
                    //----------------------------------
                    abs_b_rc = fabs(b_rc);
                    
                    if (abs_b_rc != 0.0)
                    {
                        // tau is real and will be zero if
                        // the two on-diagonal elements are
                        // equal. In this case G will be an
                        // identity matrix, and there is no
                        // point in further calculating it.
                        tau = (b_cc - b_rr) / ( 2.0 * abs_b_rc);
                        //LOGD("[JacobiSVD] tau = %1.4f",tau);

                        if (tau != 0.0)
                        {
                            Rcnt += 1;
                            //LOGD("[JacobiSVD] Rcnt = %d",Rcnt);

                            if (tau < 0 )
                            {
                                signTau = -1.0;
                            }
                            else
                            {
                                signTau = 1.0;
                            }
                            
                            theta = signTau / ( fabs(tau) + sqrt( 1.0 + pow(tau,2.0)));
                            alpha = 1.0 / sqrt(1 + pow(theta,2.0));
                            beta  = (b_rc * theta * alpha) / abs_b_rc;
                            
                            /*** Update of U and V ***/
                            // U = U * G (M by N)
                            
                            // Initialize the rotation matrix G
                            /*
                             // deprecated at 180212 by jwkang
                             G.eye(N);
                             G(r,r) = alpha;
                             G(c,c) = alpha;
                             G(r,c) = beta;
                             G(c,r) = -1.0 * beta;
                             G.showMatrix("G");
                             nextU.mult(currU, G);
                             */
                            
                            // U(:,r) = U(:,r) .* G(r,r) + U(:,c) .* G(c,r);
                            vmadd(colVectorOfU_r, colVectorOfU_c, alpha, -1.0 * beta, tempVector);
                            nextU.insert_column(tempVector, r);
                            //  U(:,c) = U(:,r) .* G(r,c) + U(:,c) .* G(c,c);
                            vmadd(colVectorOfU_r, colVectorOfU_c, beta, alpha, tempVector);
                            nextU.insert_column(tempVector, c);
                            //nextU.showMatrix("nextU");
                            
                            // The matrix V is not necessary if mIsSingularVectorMatrixExport == false
                            if (mIsSingularVectorMatrixExport)
                            {
                                G.eye(N);
                                G(r,r) = alpha;
                                G(c,c) = alpha;
                                G(r,c) = beta;
                                G(c,r) = -1.0 * beta;
                                mV.mult(currV,G); // V = V * G
                                //mV.showMatrix("mV");
                            }
                        }
                    }
                }
            }

            /*** step3: calcuate singular value from eigenvalue matrix B ***/
            Ut.transpose(nextU);
            B.mult(Ut,nextU);// B= U^T * U (N by N)
            
            triuB.triu(B);
            
            //B.showMatrix("B");
            //triuB.showMatrix("triu(B)");
            
            //    Off = sum(sum(abs(triu(B, 1))))/(N.^2);
            OffDiagCost = triuB.absSumAllElem() / pow(N,2.0);
            
            if ( (OffDiagCost / OnDiagCost) <svdIterationTol)
            {
                break;
            }
            
            if(prevOffDiagCost < OffDiagCost)
            {
                break;
            }
            prevOffDiagCost = OffDiagCost;
        }
        printf("[JacobiSVD] Spent numIter= %d for SVD finding iteration",numIter+1);
        
        B.diag(mSingularValues);
        mSingularValues.absSqrtElem();
        mNuclearNorm = mSingularValues.sumAllElem();
        //mSingularValues.showVector("mSingularValues");

        
        if (mIsSingularVectorMatrixExport)
        {
            // sigular value sorting (TBU)
            
            // matrix U,V column sorting (TBU)
            
            // Matrix U normalization by singular values
            for (unsigned int ii = 0 ; ii < N ; ii++)
            {
                double currSingularValue = mSingularValues(ii);
                if ( currSingularValue> 0.0)
                {
                    nextU.columnDivScaling(ii, currSingularValue);
                }
            }
            
            if (qropt)
            {
                // with qropt == true, M == N
                Q.sliceMatColData(tempQ.getMatrixDataInDoubleArray(),mColSize_A, N,mRowSize_A);
                //Q.showMatrix("Q"); // M by N
                //nextU.showMatrix("nextU"); // M(=N) by N
                mU.mult(Q, nextU); // M by N
            }
            else
            {
                mU = nextU;
            }
            //mU.showMatrix("U");
            //mV.showMatrix("V");
            
            diagMatrix(mSingularValues, mS);
            //mS.showMatrix("S");
        }
    }
    else
    {
        printf("[JacobiSVD] Svd is already execuated.");
    }


    /* matrix memory deallocation */
    Q.deallocate();
    tempQ.deallocate();
    R.deallocate();

    currU.deallocate();
    nextU.deallocate();
    Ut.deallocate();
    
    G.deallocate();
    currV.deallocate();
    B.deallocate();
    triuB.deallocate();
    
    colVectorOfU_c.deallocate();
    colVectorOfU_r.deallocate();
    tempVector.deallocate();
    
    endTime = clock();
    
    mSVDCompute_ms = (float)(endTime-startTime)/CLOCKS_PER_SEC*1000.0;
    printf("[JacobiSVD] SVD compute() Processing time <ms> = %3.4f",mSVDCompute_ms);

    
    return mIsSvdRun;
    
}

# ifdef __DSPLIBTEST__

void JacobiSVD::showMatrixA()
{
    mA.showMatrix("A");
}

void JacobiSVD::showMatrixS()
{
    if (mIsSvdRun && mIsSingularVectorMatrixExport)
    {
        mS.showMatrix("S");
    }
    else
    {
        LOGE("[JacobiSVD] Matrix S is not ready.");
    }
}

void JacobiSVD::showMatrixU()
{
    if (mIsSvdRun && mIsSingularVectorMatrixExport)
    {
        mU.showMatrix("U");
    }
    else
    {
        LOGE("[JacobiSVD] Matrix U is not ready.");
    }
}

void JacobiSVD::showMatrixV()
{
    if (mIsSvdRun && mIsSingularVectorMatrixExport)
    {
        mV.showMatrix("V");
    }
    else
    {
        LOGE("[JacobiSVD] Matrix V is not ready.");
    }
}

void JacobiSVD::showSingularValues()
{
    
    if (mIsSvdRun)
    {
        mSingularValues.showVector("singular values");
    }
    else
    {
        LOGD("[JacobiSVD] runSvd() is not execuated.");
    }
}

void JacobiSVD::showResult()
{
    if (mIsSvdRun)
    {
        showSingularValues();
        LOGD("[JacobiSVD] mNuclearNorm= %1.4f ",mNuclearNorm);
        if (mIsSingularVectorMatrixExport)
        {
            showMatrixU();
            showMatrixS();
            showMatrixV();
        }
    }
    else
    {
        LOGD("[JacobiSVD] SVD is not computed.");
    }

}
#endif
