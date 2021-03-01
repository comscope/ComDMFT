#ifndef INCLUDE_LINALG_LINALG_H
#define INCLUDE_LINALG_LINALG_H

#include <vector>
#include <complex>
#include <cassert>
#include <algorithm>

#include "../io/Matrix.h"
#include "../BlasLapack.h"

//TODO:    check dimensions for dgemm !!!!!

namespace linalg {
    
    template<typename Value>
    inline void mult(char transA, char transB, Value alpha, io::Matrix<Value> const& matrixA, io::Matrix<Value> const& matrixB, Value beta, io::Matrix<Value>& matrixC) {
        int M = transA == 'n' ? matrixA.I() : matrixA.J();
        int N = transB == 'n' ? matrixB.J() : matrixB.I();
        int K = transB == 'n' ? matrixB.I() : matrixB.J();
        int lda = matrixA.I();
        int ldb = matrixB.I();
        int ldc = matrixC.I();
        
        gemm(
             &transA,
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
    
    
    template<typename Value>
    inline Value trace(io::Matrix<Value> const& matrix) {
        if(matrix.I() != matrix.J())
            throw std::runtime_error("linalg::trace: matrix is not square");
        
        Value sum = .0;
        for(int i = 0; i < matrix.I(); ++i) sum += matrix(i, i);
        return sum;
    };

    template<typename Value>
    inline Value trace(io::Matrix<Value> const& matrixA, io::Matrix<Value> const& matrixB) {
        if(matrixA.J() != matrixB.I() || matrixA.I() != matrixB.J())
            throw std::runtime_error("linalg::trace: dimensions do not match");
        
        int const dimI = matrixA.I();
        int const dimJ = matrixB.J();
        
        Value temp = .0;
        for(int i = 0; i < dimI; ++i)
            for(int j = 0; j < dimJ; ++j)
                temp += matrixA(i, j)*matrixB(j, i);
        
        return temp;
    };
    
    inline void eig(char jobz, char uplo, io::rmat& matrix, io::rvec& eigen_values) {
        if(matrix.I() != matrix.J())
            throw std::runtime_error("linalg::syev: matrix is not square.");
        
        int const dim = matrix.I();
        int lwork = 3*dim - 1;
        double work[lwork]; int info;
        
        syev(
             &jobz,
             &uplo,
             &dim,
             matrix.data(),
             &dim,
             eigen_values.data(),
             work,
             &lwork,
             &info
            );
        
        assert(info == 0);
    };
    
/*
 //The QR algorithm seems to cause issues with MKL libraries for some problems. Moving to divide and conquer versions of zgesvd and zheev.
    inline void eig(char jobz, char uplo, io::cmat& matrix, io::rvec& eigen_values) {
        if(matrix.I() != matrix.J())
            throw std::runtime_error("linalg::syev: matrix is not square.");
        
        int const dim = matrix.I();
        int lwork = 2*dim - 1;
        ut::complex* work = new ut::complex[lwork];
        double rwork[3*dim - 2];
        int info;
        
        heev(
             &jobz,
             &uplo,
             &dim,
             matrix.data(),
             &dim,
             eigen_values.data(),
             work,
             &lwork,
             rwork,
             &info
            );
        
        delete[] work;
        
        assert(info == 0);
    };
*/

    inline void eig(char jobz, char uplo, io::cmat& matrix, io::rvec& eigen_values) {
        if(matrix.I() != matrix.J())
            throw std::runtime_error("linalg::syev: matrix is not square.");
        
        int const dim = matrix.I();
        int lwork = 2*dim + dim*dim;
        int lrwork = 1 + 5*dim + 2*dim*dim;
        int liwork = 3 + 5*dim;
        int info;
        
        
        ut::complex* work = new ut::complex[lwork];
        double* rwork = new double[lrwork];
        int* iwork = new int[liwork];
        
        heevd(
             &jobz,
             &uplo,
             &dim,
             matrix.data(),
             &dim,
             eigen_values.data(),
             work,
             &lwork,
             rwork,
             &lrwork,
             iwork,
             &liwork,
             &info
            );
        
        delete[] work;
        delete[] rwork;
        delete[] iwork;
        
        assert(info == 0);
    };
     

    template<typename Value>
    inline io::Matrix<Value> inv(io::Matrix<Value> const& matrix) {
        if(matrix.I() != matrix.J())
            throw std::runtime_error("linalg::inv: matrix is not square.");
        
        io::Matrix<Value> temp = matrix;
        io::Matrix<Value> invers(matrix.I(), matrix.J());
        
        for(int i = 0; i < matrix.I(); ++i)
            invers(i, i) = 1.;
        
        int ipiv[matrix.I()]; int info; int size = matrix.I();
        gesv(&size, &size, temp.data(), &size, ipiv, invers.data(), &size, &info);
        assert(info == 0);
        
        return invers;
    };
    
    
    inline double spectral_norm(io::rmat const& matrix) {
        char const jobu = 'N', jobvt = 'N';
        int const m = matrix.I(), n = matrix.J(), lda = matrix.I(), ldu = 1, ldvt = 1, lwork = std::max(1, std::max(3*std::min(m, n) + std::max(m, n), 5*std::min(m, n)));
        io::rmat a(matrix); double s[std::min(m, n)]; double work[lwork];
        int info;
        
        dgesvd_(&jobu, &jobvt, &m, &n, a.data(), &lda, s, nullptr, &ldu, nullptr, &ldvt, work, &lwork, &info);
        assert(info == 0);
        
        return s[0];
    };
    
/*
     //The QR algorithm seems to cause issues with MKL libraries for select f-shell problems. Moving to divide and conquer versions of zgesvd and zheev.
    inline double spectral_norm(io::cmat const& matrix) {
        char const jobu = 'N', jobvt = 'N';
        int const m = matrix.I(), n = matrix.J(), lda = matrix.I(), ldu = 1, ldvt = 1, lwork = std::max(1, 2*std::min(m, n) + std::max(m, n));
        io::cmat a(matrix); double s[std::min(m, n)]; std::vector<ut::complex> work(lwork); double rwork[5*std::min(m, n)];
        int info;
        
        zgesvd_(&jobu, &jobvt, &m, &n, a.data(), &lda, s, nullptr, &ldu, nullptr, &ldvt, work.data(), &lwork, rwork, &info);
        assert(info == 0);
        
        return s[0];
    };
 */


    inline double spectral_norm(io::cmat const& matrix) {
        char const jobz = 'N';
        int const m = matrix.I(), n = matrix.J(), lda = matrix.I(), ldu = 1, ldvt = 1;
        io::cmat a(matrix); std::vector<double> s(std::min(m, n));
        int mx = std::max(m,n); int mn = std::min(m,n); int lwork =  2*mn+mx;
        std::vector<ut::complex> work(lwork);
        double* rwork = new double[7*mn];
        int* iwork = new int[8*mn];
        
        int info;
        
        zgesdd_(&jobz, &m, &n, a.data(), &lda, s.data(), nullptr, &ldu, nullptr, &ldvt, work.data(), &lwork, rwork, iwork, &info);
        assert(info == 0);
        
        delete[] rwork;
        delete[] iwork;
        
        return s.at(0);
    };
 

}

#endif
