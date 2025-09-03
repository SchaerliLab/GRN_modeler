#ifndef MEX_EIG_DENSE
#define MEX_EIG_DENSE

#include "mex.h"
#include "matrix.h"
#include <eigen3/Eigen/Dense>
#include <iostream>

// A: Eigen::Map<Eigen::MatrixXd> A
// M: prhs[2]
inline Eigen::Map<Eigen::MatrixXd> mat2eig_matrix(const mxArray *M) {
    size_t nRows = mxGetM(M);
    size_t nCols = mxGetN(M);
    Eigen::Map<Eigen::MatrixXd> A(mxGetPr(M),nRows,nCols);
    return(A);
}

// A: Eigen::Map<Eigen::VectorXd>
// M: prhs[2]
inline Eigen::Map<Eigen::VectorXd> mat2eig_vector(const mxArray *M) {
    size_t nRows = mxGetM(M);
    Eigen::Map<Eigen::VectorXd> A(mxGetPr(M),nRows);
    return(A);
}

// read a double vector and convert it to size_t vector
inline Eigen::Matrix<size_t, Eigen::Dynamic, 1>
mat2eig_size_t_vector(const mxArray *M) {
    if (!mxIsDouble(M) || mxIsComplex(M)) {
        mexErrMsgIdAndTxt("MyFunc:TypeError", "Input must be a real double array.");
    }

    size_t nRows = mxGetM(M);
    Eigen::Map<const Eigen::VectorXd> temp(mxGetPr(M), static_cast<Eigen::Index>(nRows));

    Eigen::Matrix<size_t, Eigen::Dynamic, 1> result(nRows);
    for (size_t i = 0; i < nRows; ++i) {
        result(i) = static_cast<size_t>(temp(i));
    }
    return result;
}

// A: double A;
// M: prhs[2]
// str: "dt"
inline double mat2eig_double(const mxArray *M,const char *str) {
    mxArray *mxValue; 
    mxValue = mxGetField(M, 0, str);
    return((double)mxGetPr(mxValue)[0]);
}
inline int mat2eig_int(const mxArray *M,const char *str) {
    mxArray *mxValue; 
    mxValue = mxGetField(M, 0, str);
    return((int)mxGetPr(mxValue)[0]);
}
inline size_t mat2eig_size_t(const mxArray *M,const char *str) {
    mxArray *mxValue; 
    mxValue = mxGetField(M, 0, str);
    return((size_t)mxGetPr(mxValue)[0]);
}


// A: Eigen::Map<Eigen::VectorXd> A
// M: plhs[1]
inline Eigen::Map<Eigen::VectorXd> eig2mat_vector(mxArray *(&M), const size_t nRows) {
    M = mxCreateDoubleMatrix(nRows, 1, mxREAL); // Create MATLAB array
    Eigen::Map<Eigen::VectorXd> A(mxGetPr(M),nRows);
    return(A);
}

// A: Eigen::Map<Eigen::MatrixXd> A
// M: plhs[1]
inline Eigen::Map<Eigen::MatrixXd> eig2mat_matrix(mxArray *(&M), const size_t nRows, const size_t nCols) {
    M = mxCreateDoubleMatrix(nRows, nCols, mxREAL); // Create MATLAB array
    Eigen::Map<Eigen::MatrixXd> A(mxGetPr(M),nRows,nCols);
    return(A);
}

#endif // MEX_EIG_DENSE