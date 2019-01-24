#include <cstring>
#include "../../include/basen.hpp"

#include "TraceAlgebraHost.h"

#define SWITCH(arg) arg##Host
#include "../include/TraceAlgebra.h"
#undef SWITCH


using namespace TrHost;

int TrHost::Comm::iStream_ = 0;

//---------------------------------------------------------------EIGENVALUES-----------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------
TrHost::EigenValues::EigenValues(Json::Value const& jParams, Json::Value const& jBloc, double Ur0) :
dim0_(jBloc("Dimension").int64()),
dim_(jParams.is("TRUNC_DIM") ? std::min<int>(dim0_, jParams("TRUNC_DIM").int64()) : dim0_),
data_(/*Comm::alloc(dim_, 0)*/ new double[dim_]),
min_(std::numeric_limits<double>::max()) { 
    std::vector<double> eig;
    
    if(jBloc("Energies").type() == Json::Value::string_type) {   // Plutonium + Titan + json_spirit::mValue = (computer memory) gods spread cheeks to ram cock in ass !
        std::string buffer; bn::decode_b64(jBloc("Energies").string().begin(), jBloc("Energies").string().end(), back_inserter(buffer));
        eig.resize(buffer.size()/sizeof(double)); buffer.copy(reinterpret_cast<char*>(eig.data()), buffer.size());
    } else if(jBloc("Energies").type() == Json::Value::array_type) {
        auto const& jEnergies = jBloc("Energies").array(); eig.resize(jEnergies.size());
        for(std::size_t e = 0; e < jEnergies.size(); ++e) eig[e] = jEnergies[e].real64();
    } else
        throw std::runtime_error(": invalid energies value.");
    
    if(dim0_ != static_cast<int>(eig.size()))
        throw std::runtime_error(": invalid energies size.");
    
    int const N = jBloc("Filling").int64();
    double const mu = jParams("mu").real64();
    for(int i = 0; i < dim_; ++i) {
        data_[i] = eig[i] - mu*N + .5*Ur0*N*N;  //eig[i] -(mu - .5*Ur0)*N + .5*Ur0*N*(N - 1)
        min_ = std::min(min_, data_[i]);
    }
}

TrHost::EigenValues::~EigenValues() { /*Comm::free(data_, 0);*/ delete[] data_;}

	
//-------------------------------------------------------------------VECTOR------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------
TrHost::Vector::Vector(double time, EigenValues const& eig) : time_(time), exponent_(time*eig.min()), eig_(eig), data_(/*Comm::alloc(eig.dim())*/ new double[eig.dim()])  {
	for(int i = 0; i < eig.dim(); ++i) get(data_)[i] = std::exp(time*get(eig.data())[i] - exponent_);
}
		
TrHost::Vector::~Vector() { /*Comm::free(data_);*/ delete[] data_;}

	
//------------------------------------------------------------------MATRIX-------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------

TrHost::Matrix::Matrix(int size) : data_(/*Comm::alloc(size)*/ new double[size]) {}

TrHost::Matrix::Matrix(Matrix::Identity const& identity) : I_(identity.dim), J_(identity.dim), data_(/*Comm::alloc(I*J)*/ new double[I_*J_]), exponent_(.0) {
    std::memset(get(data_), 0, I_*J_*sizeof(double)); //huere memset isch das allgemein für double's ?
    for(int i = 0; i < identity.dim; ++i) get(data_)[i*(identity.dim + 1)] = 1.;
}
TrHost::Matrix::Matrix(int I, int J, std::vector<double> const& mat, int dataColMajor, int I0, int J0) : I_(I), J_(J), data_(/*Comm::alloc(I*J)*/ new double[I_*J_]), exponent_(.0) {
    if(static_cast<int>(mat.size()) != I0*J0) throw(std::runtime_error("Wrong matrix size"));
    
    for(int i = 0; i < I; ++i) for(int j = 0; j < J; ++j) data_[j + J*i] = dataColMajor ? mat[i + j*I0] : mat[j + i*J0];
}

TrHost::Matrix::Matrix(Matrix::Zero const& zero) : I_(zero.dim), J_(zero.dim), data_(/*Comm::alloc(I*J)*/ new double[I_*J_]), exponent_(.0) {
    std::memset(get(data_), 0, I_*J_*sizeof(double)); //huere memset isch das allgemein für double's ?
}
TrHost::Matrix::Matrix(double time, EigenValues const& eig) : I_(eig.dim()), J_(eig.dim()), data_(/*Comm::alloc(I*J)*/ new double[I_*J_]), exponent_(time*eig.min()) {
    std::memset(get(data_), 0, I_*J_*sizeof(double));
    for(int i = 0; i < eig.dim(); ++i) get(data_)[(eig.dim() + 1)*i] = std::exp(time*get(eig.data())[i] - exponent_);
}

TrHost::Matrix::~Matrix() { /*Comm::free(data_);*/ delete[] data_;}


void TrHost::copyEvolveL(Vector const& prop, Matrix& dest, Matrix const& source) {
    dest.I_ = source.I_; dest.J_ = source.J_; dest.exponent_ = source.exponent_ + prop.exponent(); int const inc = 1; // eigentli source.exponent_ = 0 wil basis-operator, isch aber sicherer so.
    std::memset(get(dest.data_), 0, dest.I_*dest.J_*sizeof(double));
    for(int i = 0; i < source.I_; ++i) daxpy_(&source.J_, get(prop.data()) + i, get(source.data_) + i*source.J_, &inc, get(dest.data_) + i*dest.J_, &inc);
}

	
void TrHost::mult(Matrix& dest, Matrix const& L, Matrix const& R) {
    dest.I_ = L.I_; dest.J_ = R.J_; dest.exponent_ = L.exponent_ + R.exponent_;
	char transNo = 'n'; double one = 1.; double zero = .0; 
	dgemm_(&transNo, &transNo, &R.J_, &L.I_, &L.J_, &one, get(R.data_), &R.J_, get(L.data_), &L.J_, &zero, get(dest.data_), &dest.J_);
}
	
void TrHost::evolveL(Vector const& prop, Matrix& arg) {
	arg.exponent_ += prop.exponent(); int const inc = 1;
	for(int i = 0; i < arg.I_; ++i) dscal_(&arg.J_, get(prop.data()) + i, get(arg.data_) + i*arg.J_, &inc);
}

void TrHost::trace(Zahl::Zahl* Z, Zahl::Zahl* accZ, Matrix const& matrix) {
	double sum = .0; for(int i = 0; i < matrix.I_; ++i) sum += get(matrix.data_)[(matrix.I_ + 1)*i];
    Zahl::Zahl temp(sum, matrix.exponent_); *Z = temp; *accZ += temp;
}

void TrHost::accE(double* result, Zahl::Zahl fact, Matrix const& matrix, EigenValues const& eig) {
	double sum = .0; for(int i = 0; i < matrix.I_; ++i) sum += get(eig.data())[i]*get(matrix.data_)[(matrix.I_ + 1)*i];
	*result += (fact*Zahl::Zahl(sum, matrix.exponent_)).toDouble();
}

void TrHost::norm(double* norm, Matrix const& matrix) {
	int inc = 1; int n = matrix.I_*matrix.J_;
    *norm = std::log(dnrm2_(&n, get(matrix.data_), &inc)) + matrix.exponent_;
}

void TrHost::axpy(Matrix& dest, Zahl::Zahl const& fact, Matrix const& source) {
	int const one = 1; int const N = source.I_*source.J_;
	double const x = (fact*Zahl::pow(source.exponent_)).toDouble();
	daxpy_(&N, &x, get(source.data_), &one, get(dest.data_), &one);
}

void TrHost::axpy(double* dest, double fact, Matrix const& source) {
	int const one = 1; int const N = source.I_*source.J_;
	daxpy_(&N, &fact, get(source.data_), &one, dest, &one);
}


