#ifndef LINALG
#define LINALG

#include <iostream>
#include <cstring>
#include <vector>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cerrno>
#include <utility>
#include <valarray>
#include <limits>
#include <map>
#include <complex>

#include "IO.h"


//TODO:    check dimensions for dgemm !!!!! 

extern "C" {
    void   dgemm_(const char*, const char*, const int*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
    void   zgemm_(const char*, const char*, const int*, const int*, const int*, const std::complex<double>*, const std::complex<double>*, const int*, const std::complex<double>*, const int*, const std::complex<double>*, std::complex<double>*, const int*);
    void   dsyev_(const char*, const char*, const int*, double*, const int*, double*, double*, const int*, int*);
    void   zgesv_(const int*, const int*, std::complex<double>*, const int*, int*, std::complex<double>*, const int*, int*);
};

namespace linalg {
    
    inline void gemm(char transA, char transB, double alpha, io::rmat const& matrixA, io::rmat const& matrixB, double beta, io::rmat& matrixC) {
        int M = transA == 'n' ? matrixA.I() : matrixA.J();
        int N = transB == 'n' ? matrixB.J() : matrixB.I();
        int K = transB == 'n' ? matrixB.I() : matrixB.J();
        int lda = matrixA.I();
        int ldb = matrixB.I();
        int ldc = matrixC.I();
        
        dgemm_(&transA,
               &transB,
               &M,
               &N,
               &K,
               &alpha,
               matrixA.data(),
               &lda,
               matrixB.data(),
               &ldb,
               &beta,
               matrixC.data(),
               &ldc
               );
    };
    
    inline double trace(io::rmat const& matrix) {
        if(matrix.I() != matrix.J())
            throw std::runtime_error("linalg::trace: matrix is not square");
        
        double sum = .0;
        for(int i = 0; i < matrix.I(); ++i) sum += matrix(i, i);
        return sum;
    };
    
    inline double trace(io::rmat const& matrixA, io::rmat const& matrixB) {
        if(matrixA.J() != matrixB.I() || matrixA.I() != matrixB.J())
            throw std::runtime_error("trace: dimensions do not match");
        int const dimI = matrixA.I();
        int const dimJ = matrixB.J();
        
        double temp = .0;
        for(int i = 0; i < dimI; ++i)
            for(int j = 0; j < dimJ; ++j)
                temp += matrixA(i, j)*matrixB(j, i);
        
        return temp;
    };
    
    inline void gemm(char transA, char transB, std::complex<double> alpha, io::cmat const& matrixA, io::cmat const& matrixB, std::complex<double> beta, io::cmat& matrixC) {
        int M = transA == 'n' ? matrixA.I() : matrixA.J();
        int N = transB == 'n' ? matrixB.J() : matrixB.I();
        int K = transB == 'n' ? matrixB.I() : matrixB.J();
        int lda = matrixA.I();
        int ldb = matrixB.I();
        int ldc = matrixC.I();
        
        zgemm_(&transA,
               &transB,
               &M,
               &N,
               &K,
               &alpha,
               matrixA.data(),
               &lda,
               matrixB.data(),
               &ldb,
               &beta,
               matrixC.data(),
               &ldc
               );
    };
    
    inline void syev(char jobz, char uplo, io::rmat& matrix, io::rvec& eigen_values) {
        if(matrix.I() != matrix.J())
            throw std::runtime_error("linalg::syev: matrix is not square.");
        
        int const dim = matrix.I();
        int lwork = 3*dim - 1;
        double work[lwork]; int info;
        
        dsyev_(&jobz,
               &uplo,
               &dim,
               matrix.data(),
               &dim,
               eigen_values.data(),
               work,
               &lwork,
               &info);
        
        assert(info == 0);
    };
    
    inline io::cmat inv(io::cmat const& matrix) {
        if(matrix.I() != matrix.J())
            throw std::runtime_error("linalg::inv: matrix is not square.");
        
        io::cmat temp = matrix;
        io::cmat invers;
        
        invers.resize(matrix.I(), matrix.J());
        for(int i = 0; i < matrix.I(); ++i)
            invers(i, i) = 1.;
        
        int ipiv[matrix.I()]; int info; int size = matrix.I();
        zgesv_(&size, &size, temp.data(), &size, ipiv, invers.data(), &size, &info);        
        assert(info == 0);
        
        return invers;
    };
}

#endif
