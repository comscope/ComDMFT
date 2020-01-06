#ifndef BLASLAPACK
#define BLASLAPACK

#include <complex>

#include "../ctqmc/include/Utilities.h"

extern "C" {
	double dasum_(int const*, double const*, int const*);
	double dnrm2_(int const*, double const*, int const*);
	double ddot_(int const*, double const*, int const*, double const*, int const*);
	void dswap_(int const*, double*, int const*, double*, int const*);
	void daxpy_(int const*, double const*, double const*, int const*, double*, int const*);
	void dger_(int const*, int const*, double const*, double const*, int const*, double const*, int const*, double*, int const*);
	void dscal_(int const*, double const*, double*, int const*);
	void dgemm_(char const*, char const*, int const*, int const*, int const*, double const*, double const*, int const*, double const*, int const*, double const*, double*, int const*);
	void dcopy_(int const*, double const*, int const* , double*, int const*);
	void dgemv_(char const*, int const*, int const*, double const*, double const*, int const*, double const*, int const*, double const*, double*, int const*);
	void dgesv_(int const*, int const*, double*, int const*, int*, double*, int const*, int*);
	void dgesdd_(char const*, int const*, int const*, double*, int const*, double*, double*, int const*, double*, int const*, double*, int const*, int*, int*);
    void dsyev_(char const*, char const*, int const*, double*, int const*, double*, double*, int const*, int*);

    void zcopy_(int const*, ut::complex const*, int const* , ut::complex*, int const*);
    void zgemm_(char const*, char const*, int const*, int const*, int const*, const ut::complex*, const ut::complex*, int const*, const ut::complex*, int const*, const ut::complex*, ut::complex*, int const*);
    void zgesv_(int const*, int const*, ut::complex*, int const*, int*, ut::complex*, int const*, int*);
    void zgemv_(char const*, int const*, int const*, ut::complex const*, ut::complex const*, int const*, ut::complex const*, int const*, ut::complex const*, ut::complex*, int const*);
    void zgeru_(int const*, int const*, ut::complex const*, ut::complex const*, int const*, ut::complex const*, int const*, ut::complex*, int const*);
    void zscal_(int const*, ut::complex const*, ut::complex*, int const*);
    void zswap_(int const*, ut::complex*, int const*, ut::complex*, int const*);
}

#endif
