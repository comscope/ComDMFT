#include <cstring>

#include "AlgebraHost.h"
#include "../../include/BlasLapack.h"

using namespace imp;

//---------------------------------------------------------------EIGENVALUES-----------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------
template<>
imp::Energies<AllocHost>::Energies(jsx::value const& jParams, std::vector<double> const& energies) :
dim0_(energies.size()),
dim_(jParams.is("trunc dim") ? std::min<int>(dim0_, jParams("trunc dim").int64()) : dim0_),
ln_dim_(std::log(dim_)),
data_(dim_),
min_(std::numeric_limits<double>::max()) {
    for(int i = 0; i < dim_; ++i) {
        min_ = std::min(min_, energies[i]);
        data_.ptr()[i] = energies[i];
    }
}

template<>
imp::Energies<AllocHost>::~Energies() {
}

	
//-------------------------------------------------------------------VECTOR------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------
template<>
imp::Vector<AllocHost>::Vector(double time, Energies<AllocHost> const& energies) : time_(time), exponent_(time*energies.min()), energies_(energies), data_(energies.dim())  {
	for(int i = 0; i < energies.dim(); ++i) data_.ptr()[i] = std::exp(time*energies.data().ptr()[i] - exponent_);
}

template<>
imp::Vector<AllocHost>::~Vector() {
}

	
//------------------------------------------------------------------MATRIX-------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------
template<>
imp::Matrix<AllocHost>::Matrix(int size) : data_(size) {
}
template<>
imp::Matrix<AllocHost>::Matrix(Matrix<AllocHost>::Identity const& identity) : I_(identity.dim), J_(identity.dim), data_(I_*J_), exponent_(.0) {
    std::memset(data_.ptr(), 0, I_*J_*sizeof(double)); //huere memset isch das allgemein für double's ?
    for(int i = 0; i < identity.dim; ++i) data_.ptr()[i*(identity.dim + 1)] = 1.;
}
template<>
imp::Matrix<AllocHost>::Matrix(Matrix<AllocHost>::Zero const& zero) : I_(zero.dim), J_(zero.dim), data_(I_*J_), exponent_(.0) {
    std::memset(data_.ptr(), 0, I_*J_*sizeof(double)); //huere memset isch das allgemein für double's ?
}
template<>
imp::Matrix<AllocHost>::Matrix(int I, int J, io::rmat const& matrix) : I_(I), J_(J), data_(I_*J_), exponent_(.0) {
    for(int i = 0; i < I; ++i)
        for(int j = 0; j < J; ++j)
            data_.ptr()[j + J*i] = matrix(i, j);
}

template<>
imp::Matrix<AllocHost>::~Matrix() {
}


void imp::copyEvolveL(Matrix<AllocHost>& dest, Vector<AllocHost> const& prop, Matrix<AllocHost> const& source, Batcher<AllocHost>& batcher) {
    dest.I() = source.I(); dest.J() = source.J(); dest.exponent() = source.exponent() + prop.exponent(); int const inc = 1; // eigentli source.exponent_ = 0 wil basis-operator, isch aber sicherer so.
    std::memset(dest.data().ptr(), 0, dest.I()*dest.J()*sizeof(double));
    for(int i = 0; i < source.I(); ++i) daxpy_(&source.J(), prop.data().ptr() + i, source.data().ptr() + i*source.J(), &inc, dest.data().ptr() + i*dest.J(), &inc);
}

void imp::mult(Matrix<AllocHost>& dest, Matrix<AllocHost> const& L, Matrix<AllocHost> const& R, Batcher<AllocHost>& batcher) {
    dest.I() = L.I(); dest.J() = R.J(); dest.exponent() = L.exponent() + R.exponent();
	char transNo = 'n'; double one = 1.; double zero = .0; 
	dgemm_(&transNo, &transNo, &R.J(), &L.I(), &L.J(), &one, R.data().ptr(), &R.J(), L.data().ptr(), &L.J(), &zero, dest.data().ptr(), &dest.J());
}
	
void imp::evolveL(Vector<AllocHost> const& prop, Matrix<AllocHost>& arg, Batcher<AllocHost>& batcher) {
	arg.exponent() += prop.exponent(); int const inc = 1;
	for(int i = 0; i < arg.I(); ++i) dscal_(&arg.J(), prop.data().ptr() + i, arg.data().ptr() + i*arg.J(), &inc);
}

void imp::trace(ut::Zahl* Z, ut::Zahl* accZ, Matrix<AllocHost> const& matrix, Batcher<AllocHost>& batcher) {
	double sum = .0; for(int i = 0; i < matrix.I(); ++i) sum += matrix.data().ptr()[(matrix.I() + 1)*i];
    ut::Zahl temp(sum, matrix.exponent()); if(Z) *Z = temp; if(accZ) *accZ += temp;
}

void imp::traceAtB(ut::Zahl* Z, ut::Zahl* accZ, Matrix<AllocHost> const& At, Matrix<AllocHost> const& B, Batcher<AllocHost>& batcher) {
    if(At.I() != B.I() || At.J() != B.J()) throw std::runtime_error("imp::traceAtB: scheisse");
    int const n = At.I()*At.J(); int const inc = 1;
    double sum = ddot_(&n, At.data().ptr(), &inc, B.data().ptr(), &inc);
    ut::Zahl temp(sum, At.exponent() + B.exponent()); if(Z) *Z = temp; if(accZ) *accZ += temp;
}

void imp::norm(double* norm, Matrix<AllocHost> const& matrix, Batcher<AllocHost>& batcher) {
	int inc = 1; int n = matrix.I()*matrix.J();
    *norm = std::log(dnrm2_(&n, matrix.data().ptr(), &inc)) + matrix.exponent();
}

void imp::density_matrix(Matrix<AllocHost>& dest, ut::Zahl const& fact, Matrix<AllocHost> const& Bt, Vector<AllocHost> const& prop, Matrix<AllocHost> const& A, Matrix<AllocHost>& buffer, Batcher<AllocHost>& batcher) {
    buffer.I() = A.I(); buffer.J() = Bt.I(); buffer.exponent() = A.exponent() + Bt.exponent();
    char transNo = 'n'; char transYes = 't'; double one = 1.; double zero = .0;
    dgemm_(&transYes, &transNo, &Bt.I(), &A.I(), &A.J(), &one, Bt.data().ptr(), &Bt.J(), A.data().ptr(), &A.J(), &zero, buffer.data().ptr(), &buffer.J());
    
    buffer.exponent() += prop.exponent(); double const deltaTime = -prop.time(); // this is confusing, change time -> -time
    auto d = buffer.data().ptr(); auto const p = prop.data().ptr(); auto const e = prop.energies().data().ptr();
    for(int i = 0; i < buffer.I(); ++i)
        for(int j = 0; j < buffer.J(); ++j) {
            double const deltaE = e[j] - e[i]; double const delta = deltaTime*deltaE;
            d[j + buffer.J()*i] *= std::abs(delta) > 1.e-7 ? (p[i] - p[j])/deltaE : deltaTime/2.*(p[i] + p[j] + delta/2.*(p[j] - p[i])); //approximation is symmetric
        }
    
    dest.I() = buffer.I(); dest.J() = buffer.J();
    axpy(dest, fact, buffer, batcher);
}

void imp::axpy(Matrix<AllocHost>& dest, ut::Zahl const& fact, Matrix<AllocHost> const& source, Batcher<AllocHost>& batcher) {
	int const one = 1; int const N = source.I()*source.J();
	double const x = (fact*ut::exp(source.exponent())).to_double();
	daxpy_(&N, &x, source.data().ptr(), &one, dest.data().ptr(), &one);
}

void imp::axpy(double* dest, double fact, Matrix<AllocHost> const& source) {
	int const one = 1; int const N = source.I()*source.J();
	daxpy_(&N, &fact, source.data().ptr(), &one, dest, &one);
}



