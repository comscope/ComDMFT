#ifndef BATH_ALGEBRA_H
#define BATH_ALGEBRA_H

#include <vector>

#include "../Utilities.h"


namespace bath {
    
    template<typename HybVal>
    struct Matrix {
        Matrix() = delete;
        explicit Matrix(int dim) : dim_(dim), data_(dim_*dim_) {};
        Matrix(Matrix const&) = delete;
        Matrix(Matrix&&) = default;
        Matrix& operator=(Matrix const&) = delete;
        Matrix& operator=(Matrix&&) = default;
        ~Matrix() = default;
        
        int dim() const { return dim_;};
        HybVal& at(int i, int j) { return data_[i + dim_*j];};
        HybVal const& at(int i, int j) const { return data_[i + dim_*j];};
        HybVal* data() { return data_.data();};
        HybVal const* data() const { return data_.data();};
        HybVal* data(int i, int j) { return data_.data() + i + j*dim_;};
        HybVal const* data(int i, int j) const { return data_.data() + i + j*dim_;};
        
    private:
        int dim_;
        std::vector<HybVal> data_;
    };

    
    double xTy(int const* n, double const* dx, int const* incx, double const* dy, int const* incy) { return ddot_(n, dx, incx, dy, incy); }
    void copy(int const* n, double const* dx, int const* incx, double* dy, int const* incy) { dcopy_(n, dx, incx, dy, incy);};
    void gesv(int const* n, int const* nrhs, double* a, int const* lda, int* ipiv, double* b, int const* ldb, int* info) { dgesv_(n, nrhs, a, lda, ipiv, b, ldb, info); }
    void gemv(const char* trans, int const* m, int const* n, double const* alpha, double const* a, int const* lda, double const* x, int const* incx, double const* beta, double* y, int const* incy) { dgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy); }
    void axyT(int const* m, int const* n, double const* alpha, double const* x, int const* incx, double const* y, int const* incy, double* a, int const* lda) { dger_(m, n, alpha, x, incx, y, incy, a, lda); }
    void scal(int const* n, double const* da, double* dx, int const* incx) { dscal_(n, da, dx, incx); }
    void swap(int const* n, double* dx, int const* incx, double* dy, int const* incy) { dswap_(n, dx, incx, dy, incy); }
    
    ut::complex xTy(int const* n, ut::complex const* dx, int const* incx, ut::complex const* dy, int const* incy) {
        int const N = *n, ix = *incx, iy = *incy; ut::complex temp = .0;
        
        if(ix == 1 && iy == 1)
            for(int i = 0; i < N; ++i) temp += dx[i]*dy[i];
        else
            for(int i = 0; i < N; ++i) temp += dx[i*ix]*dy[i*iy];
        
        return temp;
    }
    void copy(int const* n, ut::complex const* dx, int const* incx, ut::complex* dy, int const* incy) { zcopy_(n, dx, incx, dy, incy);};
    void gesv(int const* n, int const* nrhs, ut::complex* a, int const* lda, int* ipiv, ut::complex* b, int const* ldb, int* info) { zgesv_(n, nrhs, a, lda, ipiv, b, ldb, info); }
    void gemv(const char* trans, int const* m, int const* n, ut::complex const* alpha, ut::complex const* a, int const* lda, ut::complex const* x, int const* incx, ut::complex const* beta, ut::complex* y, int const* incy) { zgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy); }
    void axyT(int const* m, int const* n, ut::complex const* alpha, ut::complex const* x, int const* incx, ut::complex const* y, int const* incy, ut::complex* a, int const* lda) { zgeru_(m, n, alpha, x, incx, y, incy, a, lda); }
    void scal(int const* n, ut::complex const* da, ut::complex* dx, int const* incx) { zscal_(n, da, dx, incx); }
    void swap(int const* n, ut::complex* dx, int const* incx, ut::complex* dy, int const* incy) { zswap_(n, dx, incx, dy, incy); }
    
}


#endif
