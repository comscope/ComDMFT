#ifndef INCLUDE_BLASLAPACK_H
#define INCLUDE_BLASLAPACK_H

#include <complex>

#include "../ctqmc/include/Utilities.h"

//--------------------------------------------------------------------- real --------------------------------------------------------------------

extern "C" {
	double dnrm2_(int const*, double const*, int const*);
	double ddot_(int const*, double const*, int const*, double const*, int const*);
	void dswap_(int const*, double*, int const*, double*, int const*);
    void dcopy_(int const*, double const*, int const* , double*, int const*);
    void dscal_(int const*, double const*, double*, int const*);
	void daxpy_(int const*, double const*, double const*, int const*, double*, int const*);
	void dger_(int const*, int const*, double const*, double const*, int const*, double const*, int const*, double*, int const*);
    void dgemv_(char const*, int const*, int const*, double const*, double const*, int const*, double const*, int const*, double const*, double*, int const*);
    void dgemm_(char const*, char const*, int const*, int const*, int const*, double const*, double const*, int const*, double const*, int const*, double const*, double*, int const*);
    void dsyev_(char const*, char const*, int const*, double*, int const*, double*, double*, int const*, int*);
    void dgesv_(int const*, int const*, double*, int const*, int*, double*, int const*, int*);
    void dgesdd_(char const*, int const*, int const*, double*, int const*, double*, double*, int const*, double*, int const*, double*, int const*, int*, int*);
    void dgesvd_(char const*, char const*, int const*, int const*, double*, int const*, double*, double*, int const*, double*, int const*, double*, int const*, int*);
}

inline double nrm2(int const* n, double const* dx, int const* inc) {
    return dnrm2_(n, dx, inc);
}
inline double dotu(int const* n, double const* dx, int const* incx, double const* dy, int const* incy) {
    return ddot_(n, dx, incx, dy, incy);
}
inline double dotc(int const* n, double const* dx, int const* incx, double const* dy, int const* incy) {
    return ddot_(n, dx, incx, dy, incy);
}
inline void swap(int const* n, double* dx, int const* incx, double* dy, int const* incy) {
    dswap_(n, dx, incx, dy, incy);
}
inline void copy(int const* n, double const* dx, int const* incx, double* dy, int const* incy) {
    dcopy_(n, dx, incx, dy, incy);
}
inline void scal(int const* n, double const* da, double* dx, int const* incx) {
    dscal_(n, da, dx, incx);
}
inline void axpy(int const* n, double const* a, double const* dx, int const* incx, double* dy, int const* incy) {
    daxpy_(n, a, dx, incx, dy, incy);
}
inline void geru(int const* m, int const* n, double const* alpha, double const* x, int const* incx, double const* y, int const* incy, double* a, int const* lda) {
    dger_(m, n, alpha, x, incx, y, incy, a, lda);
}
inline void gemv(const char* trans, int const* m, int const* n, double const* alpha, double const* a, int const* lda, double const* x, int const* incx, double const* beta, double* y, int const* incy) {
    dgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}
inline void gemm(char const* transa, char const* transb, int const* m, int const* n, int const* k, double const* alpha, double const* a, int const* lda, double const* b, int const* ldb, double const* beta, double* c, int const* ldc) {
    dgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}
inline void syev(char const* jobz, char const* uplo, int const* n, double* a, int const* lda, double* w, double* work, int const* lwork, int* info) {
    dsyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
}
inline void gesv(int const* n, int const* nrhs, double* a, int const* lda, int* ipiv, double* b, int const* ldb, int* info) {
    dgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
}
inline void gesvd(char const* jobu, char const* jobvt, int const* m, int const* n, double* a, int const* lda, double* s, double* u, int const* ldu, double* vt, int const* ldvt, double* work, int const* lwork, int* info) {
    dgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
}


//--------------------------------------------------------------------- complex --------------------------------------------------------------------

extern "C" {
    double dznrm2_(int const*, ut::complex const*, int const*);
    //ut::complex zdotc_(int const*, ut::complex const*, int const*, ut::complex const*, int const*);
    void zswap_(int const*, ut::complex*, int const*, ut::complex*, int const*);
    void zcopy_(int const*, ut::complex const*, int const* , ut::complex*, int const*);
    void zscal_(int const*, ut::complex const*, ut::complex*, int const*);
    void zaxpy_(int const*, ut::complex const*, ut::complex const*, int const*, ut::complex*, int const*);
    void zgeru_(int const*, int const*, ut::complex const*, ut::complex const*, int const*, ut::complex const*, int const*, ut::complex*, int const*);
    void zgemv_(char const*, int const*, int const*, ut::complex const*, ut::complex const*, int const*, ut::complex const*, int const*, ut::complex const*, ut::complex*, int const*);
    void zgemm_(char const*, char const*, int const*, int const*, int const*, const ut::complex*, const ut::complex*, int const*, const ut::complex*, int const*, const ut::complex*, ut::complex*, int const*);
    void zheev_(char const*, char const*, int const*, ut::complex*, int const* lda, double*, ut::complex*, int*, double*, int*);
    void zheevd_(char const*, char const*, int const*, ut::complex*, int const*, double*, ut::complex*, int*, double*, int*, int*, int*, int*);
    void zgesv_(int const*, int const*, ut::complex*, int const*, int*, ut::complex*, int const*, int*);
    void zgesvd_(char const*, char const*, int const*, int const*, ut::complex*, int const*, double*, ut::complex*, int const*, ut::complex*, int const*, ut::complex*, int const*, double*, int*);
    void zgesdd_(char const*, int const*, int const*, ut::complex*, int const*, double*, ut::complex*, int const*, ut::complex*, int const*, ut::complex*, int const*, double*, int*, int*);
}


inline double nrm2(int const* n, ut::complex const* x, int const* inc) {
    return dznrm2_(n, x, inc);
}
inline ut::complex dotu(int const* n, ut::complex const* dx, int const* incx, ut::complex const* dy, int const* incy) {
    int const N = *n, ix = *incx, iy = *incy; ut::complex temp = .0;
    
    if(ix == 1 && iy == 1)
        for(int i = 0; i < N; ++i) temp += dx[i]*dy[i];
    else
        for(int i = 0; i < N; ++i) temp += dx[i*ix]*dy[i*iy];
    
    return temp;
}
inline ut::complex dotc(int const* n, ut::complex const* dx, int const* incx, ut::complex const* dy, int const* incy) {
    int const N = *n, ix = *incx, iy = *incy; ut::complex temp = .0;
    
    if(ix == 1 && iy == 1)
        for(int i = 0; i < N; ++i) temp += std::conj(dx[i])*dy[i];
    else
        for(int i = 0; i < N; ++i) temp += std::conj(dx[i*ix])*dy[i*iy];
    
    return temp;
}
inline void swap(int const* n, ut::complex* dx, int const* incx, ut::complex* dy, int const* incy) {
    zswap_(n, dx, incx, dy, incy);
}
inline void copy(int const* n, ut::complex const* dx, int const* incx, ut::complex* dy, int const* incy) {
    zcopy_(n, dx, incx, dy, incy);
}
inline void scal(int const* n, double const* alpha, ut::complex* x, int const* inc) {
    if(*inc != 1) throw std::runtime_error("scal: incx is not 1");
    int const n2 = 2**n; dscal_(&n2, alpha, reinterpret_cast<double*>(x), inc);
}
inline void scal(int const* n, ut::complex const* alpha, ut::complex* x, int const* inc) {
    zscal_(n, alpha, x, inc);
}
inline void axpy(int const* n, double const* a, ut::complex const* dx, int const* incx, ut::complex* dy, int const* incy) {
    if(*incx != 1 || *incy != 1) throw std::runtime_error("axpy: incx or incy is not 1");
    int n2 = 2**n; daxpy_(&n2, a, reinterpret_cast<double const*>(dx), incx, reinterpret_cast<double*>(dy), incy);
}
inline void axpy(int const* n, ut::complex const* a, ut::complex const* dx, int const* incx, ut::complex* dy, int const* incy) {
    zaxpy_(n, a, dx, incx, dy, incy);
}
inline void geru(int const* m, int const* n, ut::complex const* alpha, ut::complex const* x, int const* incx, ut::complex const* y, int const* incy, ut::complex* a, int const* lda) {
    zgeru_(m, n, alpha, x, incx, y, incy, a, lda);
}
inline void gemv(const char* trans, int const* m, int const* n, ut::complex const* alpha, ut::complex const* a, int const* lda, ut::complex const* x, int const* incx, ut::complex const* beta, ut::complex* y, int const* incy) {
    zgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}
inline void gemm(char const* transa, char const* transb, int const* m, int const* n, int const* k, ut::complex const* alpha, ut::complex const* a, int const* lda, ut::complex const* b, int const* ldb, ut::complex const* beta, ut::complex* c, int const* ldc) {
    zgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void heev(char const* jobz, char const* uplo, int const* n, ut::complex* a, int const* lda, double* w, ut::complex* work, int* lwork, double* rwork, int* info) {
    zheev_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
}
inline void heevd(char const* jobz, char const* uplo, int const* n, ut::complex* a, int const* lda, double* w, ut::complex* work, int* lwork, double* rwork, int* lrwork, int* iwork, int* liwork, int* info) {
    zheevd_(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork,iwork, liwork, info);
}
inline void gesv(int const* n, int const* nrhs, ut::complex* a, int const* lda, int* ipiv, ut::complex* b, int const* ldb, int* info) {
    zgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
}
inline void gesvd(char const* jobu, char const* jobvt, int const* m, int const* n, ut::complex* a, int const* lda, double* s, ut::complex* u, int const* ldu, ut::complex* vt, int const* ldvt, ut::complex* work, int const* lwork, double* rwork, int* info) {
    zgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
}






#endif
