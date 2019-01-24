#include <cstring>

#include "TraceAlgebraHost.h"
#include "../include/TraceAlgebra.h"


using namespace tr;

//std::int64_t tr::counter;

int tr::Comm::iStream_ = 0;

//---------------------------------------------------------------EIGENVALUES-----------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------
tr::Energies::Energies(jsx::value const& jParams, std::vector<double> const& eig) :
dim0_(eig.size()),
dim_(jParams.is("trunc dim") ? std::min<int>(dim0_, jParams("trunc dim").int64()) : dim0_),
data_(new double[dim_]),
min_(std::numeric_limits<double>::max()) {
    for(int i = 0; i < dim_; ++i) {
        min_ = std::min(min_, eig[i]);
        data_[i] = eig[i];
    }
    
    //counter += dim_;
}

tr::Energies::~Energies() {
    delete[] data_;
    
    //counter -= dim_;
}

	
//-------------------------------------------------------------------VECTOR------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------
tr::Vector::Vector(double time, Energies const& eig) : /*dim_(eig.dim()),*/ time_(time), exponent_(time*eig.min()), /*eig_(eig),*/ data_(new double[eig.dim()])  {
	for(int i = 0; i < eig.dim(); ++i) get(data_)[i] = std::exp(time*get(eig.data())[i] - exponent_);
    
    //counter += dim_;
}
		
tr::Vector::~Vector() {
    delete[] data_;

    //counter -= dim_;
}

	
//------------------------------------------------------------------MATRIX-------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------

tr::Matrix::Matrix(int size) : /*dim_(size),*/ data_(new double[size]) {
    //counter += dim_;
}

tr::Matrix::Matrix(Matrix::Identity const& identity) : /*dim_(identity.dim*identity.dim),*/ I_(identity.dim), J_(identity.dim), data_(new double[I_*J_]), exponent_(.0) {
    std::memset(get(data_), 0, I_*J_*sizeof(double)); //huere memset isch das allgemein für double's ?
    for(int i = 0; i < identity.dim; ++i) get(data_)[i*(identity.dim + 1)] = 1.;
    
    //counter += dim_;
}
tr::Matrix::Matrix(int I, int J, io::rmat const& matrix) : /*dim_(I*J),*/ I_(I), J_(J), data_(new double[I_*J_]), exponent_(.0) {
    for(int i = 0; i < I; ++i)
        for(int j = 0; j < J; ++j)
            data_[j + J*i] = matrix(i, j);
    
    //counter += dim_;
}

tr::Matrix::Matrix(Matrix::Zero const& zero) : /*dim_(zero.dim*zero.dim),*/ I_(zero.dim), J_(zero.dim), data_(new double[I_*J_]), exponent_(.0) {
    std::memset(get(data_), 0, I_*J_*sizeof(double)); //huere memset isch das allgemein für double's ?
    
    //counter += dim_;
}
tr::Matrix::Matrix(double time, Energies const& eig) : /*dim_(eig.dim()*eig.dim()),*/ I_(eig.dim()), J_(eig.dim()), data_(new double[I_*J_]), exponent_(time*eig.min()) {
    std::memset(get(data_), 0, I_*J_*sizeof(double));
    for(int i = 0; i < eig.dim(); ++i) get(data_)[(eig.dim() + 1)*i] = std::exp(time*get(eig.data())[i] - exponent_);
    
    //counter += dim_;
}

tr::Matrix::~Matrix() {
    delete[] data_;

    //counter -= dim_;
}


void tr::copyEvolveL(Vector const& prop, Matrix& dest, Matrix const& source) {
    dest.I_ = source.I_; dest.J_ = source.J_; dest.exponent_ = source.exponent_ + prop.exponent(); int const inc = 1; // eigentli source.exponent_ = 0 wil basis-operator, isch aber sicherer so.
    std::memset(get(dest.data_), 0, dest.I_*dest.J_*sizeof(double));
    for(int i = 0; i < source.I_; ++i) daxpy_(&source.J_, get(prop.data()) + i, get(source.data_) + i*source.J_, &inc, get(dest.data_) + i*dest.J_, &inc);
}

	
void tr::mult(Matrix& dest, Matrix const& L, Matrix const& R) {
    dest.I_ = L.I_; dest.J_ = R.J_; dest.exponent_ = L.exponent_ + R.exponent_;
	char transNo = 'n'; double one = 1.; double zero = .0; 
	dgemm_(&transNo, &transNo, &R.J_, &L.I_, &L.J_, &one, get(R.data_), &R.J_, get(L.data_), &L.J_, &zero, get(dest.data_), &dest.J_);
}
	
void tr::evolveL(Vector const& prop, Matrix& arg) {
	arg.exponent_ += prop.exponent(); int const inc = 1;
	for(int i = 0; i < arg.I_; ++i) dscal_(&arg.J_, get(prop.data()) + i, get(arg.data_) + i*arg.J_, &inc);
}

void tr::trace(za::Zahl* Z, za::Zahl* accZ, Matrix const& matrix) {
	double sum = .0; for(int i = 0; i < matrix.I_; ++i) sum += get(matrix.data_)[(matrix.I_ + 1)*i];
    za::Zahl temp(sum, matrix.exponent_); *Z = temp; *accZ += temp;
}

void tr::accE(double* result, za::Zahl fact, Matrix const& matrix, Energies const& eig) {
	double sum = .0; for(int i = 0; i < matrix.I_; ++i) sum += get(eig.data())[i]*get(matrix.data_)[(matrix.I_ + 1)*i];
	*result += (fact*za::Zahl(sum, matrix.exponent_)).toDouble();
}

void tr::norm(double* norm, Matrix const& matrix) {
	int inc = 1; int n = matrix.I_*matrix.J_;
    *norm = std::log(dnrm2_(&n, get(matrix.data_), &inc)) + matrix.exponent_;
}

void tr::axpy(Matrix& dest, za::Zahl const& fact, Matrix const& source) {
	int const one = 1; int const N = source.I_*source.J_;
	double const x = (fact*za::pow(source.exponent_)).toDouble();
	daxpy_(&N, &x, get(source.data_), &one, get(dest.data_), &one);
}

void tr::axpy(double* dest, double fact, Matrix const& source) {
	int const one = 1; int const N = source.I_*source.J_;
	daxpy_(&N, &fact, get(source.data_), &one, dest, &one);
}


