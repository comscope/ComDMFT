#ifndef BLASLAPACK
#define BLASLAPACK

#include <complex>

extern "C" {
	double dasum_(int const*, double const*, int const*);
	double dnrm2_(int const*, double const*, int const*);
	double ddot_(const int*, const double*, const int*, const double*, const int*);
	void   dswap_(const int*, double*, int const*, double*, int const*);
	void   daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
	void   dger_(const int*, const int*, double const*, double const*, int const*, double const*, int const*, double*, int const*);
	void   dscal_(int const*, double const*, double*, int const*);
	void   dgemm_(const char*, const char*, int const*, int const*, int const*, double const*, double const*, int const*, double const*, int const*, double const*, double*, int const*);
	void   dcopy_(int const*, double const*, int const* , double*, int const*);
	void   dgemv_(const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
	void   dgesv_(const int*, const int*, double*, const int*, int*, double*, const int*, int*);
	void   dgesdd_(const char*, const int*, const int*, double*, int const*, double*, double*, int const*, double*, int const*, double*, int const*, int*, int*);

    void   zgemm_(const char*, const char*, const int*, const int*, const int*, const std::complex<double>*, const std::complex<double>*, const int*, const std::complex<double>*, const int*, const std::complex<double>*, std::complex<double>*, const int*);
}

#endif
