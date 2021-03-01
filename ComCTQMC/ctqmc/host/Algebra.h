#ifndef CTQMC_HOST_ALGEBRA_H
#define CTQMC_HOST_ALGEBRA_H

#include <stdexcept>
#include <vector>

#include "../include/Utilities.h"
#include "../include/impurity/Algebra.h"

#include "../../include/BlasLapack.h"
#include "../../include/JsonX.h"
#include "../../include/io/Matrix.h"

namespace imp {

    struct Host {};    
    
    template<typename Value>
    struct Batcher<Host, Value> : itf::Batcher<Value> {
        Batcher() = delete;
        Batcher(std::size_t) {};
        Batcher(Batcher const&) = delete;
        Batcher(Batcher&&) = delete;
        Batcher& operator=(Batcher const&) = delete;
        Batcher& operator=(Batcher&&) = delete;
        ~Batcher() = default;
        
        int is_ready() { return 1;};
        void launch() {};
    };
    
    
    template<>
    struct Energies<Host> {
        Energies() = delete;
        Energies(jsx::value const& jParams, std::vector<double> const& energies) :
        dim0_(energies.size()),
        dim_(jParams.is("trunc dim") ? std::min<int>(dim0_, jParams("trunc dim").int64()) : dim0_),
        ln_dim_(std::log(dim_)),
        data_(new double[dim_]),
        min_(std::numeric_limits<double>::max()) {
            for(int i = 0; i < dim_; ++i) {
                min_ = std::min(min_, energies[i]);
                data_[i] = energies[i];
            }
        };
        Energies(Energies const&) = delete;
        Energies(Energies&&) = delete;
        Energies& operator=(Energies const&) = delete;
        Energies& operator=(Energies&&) = delete;
        ~Energies() {
            delete[] data_;
        }
        
        int const& dim0() const { return dim0_;};
        int const& dim() const { return dim_;}
        double const& ln_dim() const { return ln_dim_;};
        double* data() { return data_;}
        double const* data() const { return data_;}
        double const& min() const { return min_;}
        
    private:
        int const dim0_;
        int const dim_;
        double const ln_dim_;
        double* data_;
        double min_;
    };
    
    template<>
    struct Vector<Host> {
        Vector() = delete;
        Vector(double time, Energies<Host> const& energies) :
        time_(time),
        exponent_(time*energies.min()),
        data_(new double[energies.dim()])  {
            for(int i = 0; i < energies.dim(); ++i) data_[i] = std::exp(time*energies.data()[i] - exponent_);
        };
        Vector(Vector const&) = delete;
        Vector(Vector&&) = delete;
        Vector& operator=(Vector const&) = delete;
        Vector& operator=(Vector&&) = delete;
        ~Vector() {
            delete[] data_;
        };
        
        double const& time() const { return time_;};
        double const& exponent() const { return exponent_;}
        double* data() { return data_;}
        double const* data() const { return data_;}
        
    private:
        double const time_;
        double const exponent_;
        double* data_;
    };
    
    template<typename Value>
    struct Matrix<Host, Value> {
        struct Identity { Identity(int d) : dim(d) {}; int const dim;};
        struct Zero { Zero(int d) : dim(d) {}; int const dim;};
        
        Matrix() = delete;
        Matrix(int size):
        data_(new Value[size]) {
        };
        Matrix(Identity const& identity) :
        I_(identity.dim), J_(identity.dim),
        data_(new Value[I_*J_]),
        exponent_(.0) {
            std::memset(data_, 0, I_*J_*sizeof(Value)); //huere memset isch das allgemein für double's ?
            for(int i = 0; i < identity.dim; ++i) data_[i*(identity.dim + 1)] = 1.;
        };
        Matrix(Zero const& zero) :
        I_(zero.dim), J_(zero.dim),
        data_(new Value[I_*J_]),
        exponent_(.0) {
            std::memset(data_, 0, I_*J_*sizeof(Value)); //huere memset isch das allgemein für double's ?
        };
        Matrix(int I, int J, io::Matrix<Value> const& matrix) :
        I_(I), J_(J),
        data_(new Value[I_*J_]),
        exponent_(.0) {
            for(int i = 0; i < I; ++i)
                for(int j = 0; j < J; ++j)
                    data_[j + J*i] = matrix(i, j);
        };
        Matrix(Matrix const&) = delete;
        Matrix(Matrix&&) = delete;
        Matrix& operator=(Matrix const&) = delete;
        Matrix& operator=(Matrix&&) = delete;
        ~Matrix() {
            delete[] data_;
        }
        
        int& I() { return I_;}
        int& J() { return J_;}
        int const& I() const { return I_;}
        int const& J() const { return J_;}
        Value* data() { return data_;}
        Value const* data() const { return data_;}
        double& exponent() { return exponent_;}
        double const& exponent() const { return exponent_;}
        
    private:
        int I_, J_;
        Value* data_;
        double exponent_;
    };
    
    
    template<typename Value>
    void copyEvolveL(Matrix<Host, Value>& dest, Vector<Host> const& prop, Matrix<Host, Value> const& source, itf::Batcher<Value>& batcher) {
        dest.I() = source.I(); dest.J() = source.J(); dest.exponent() = source.exponent() + prop.exponent(); int const inc = 1; // eigentli source.exponent_ = 0 wil basis-operator, isch aber sicherer so.
        std::memset(dest.data(), 0, dest.I()*dest.J()*sizeof(Value));
        for(int i = 0; i < source.I(); ++i) axpy(&source.J(), prop.data() + i, source.data() + i*source.J(), &inc, dest.data() + i*dest.J(), &inc);
    };
    
    template<typename Value>
    void mult(Matrix<Host, Value>& dest, Matrix<Host, Value> const& L, Matrix<Host, Value> const& R, itf::Batcher<Value>& batcher) {
        dest.I() = L.I(); dest.J() = R.J(); dest.exponent() = L.exponent() + R.exponent();
        char transNo = 'n'; Value one = 1.; Value zero = .0;
        gemm(&transNo, &transNo, &R.J(), &L.I(), &L.J(), &one, R.data(), &R.J(), L.data(), &L.J(), &zero, dest.data(), &dest.J());
    };
    
    template<typename Value>
    void evolveL(Vector<Host> const& prop, Matrix<Host, Value>& arg, itf::Batcher<Value>& batcher) {
        arg.exponent() += prop.exponent(); int const inc = 1;
        for(int i = 0; i < arg.I(); ++i) scal(&arg.J(), prop.data() + i, arg.data() + i*arg.J(), &inc);
    };
    
    template<typename Value>
    void trace(ut::Zahl<Value>* Z, ut::Zahl<Value>* accZ, Matrix<Host, Value> const& matrix, itf::Batcher<Value>& batcher) {
        Value sum = .0; for(int i = 0; i < matrix.I(); ++i) sum += matrix.data()[(matrix.I() + 1)*i];
        ut::Zahl<Value> temp(sum, matrix.exponent()); if(Z) *Z = temp; if(accZ) *accZ += temp;
    };
    
    template<typename Value>
    void traceAtB(ut::Zahl<Value>* Z, ut::Zahl<Value>* accZ, Matrix<Host, Value> const& A, Matrix<Host, Value> const& B, itf::Batcher<Value>& batcher) {
        if(A.I() != B.I() || A.J() != B.J()) throw std::runtime_error("traceAtB: scheisse");
        int const n = A.I()*A.J(); int const inc = 1;
        Value sum = dotc(&n, A.data(), &inc, B.data(), &inc);  // trace(AtB) = trace(BAt) and A, B are row major => trace(BAt) = < A.data(), B.data() >
        ut::Zahl<Value> temp(sum, A.exponent() + B.exponent()); if(Z) *Z = temp; if(accZ) *accZ += temp;
    };
    
    template<typename Value>
    void norm(double* norm, Matrix<Host, Value> const& matrix, itf::Batcher<Value>& batcher) {
        int inc = 1; int n = matrix.I()*matrix.J();
        *norm = std::log(nrm2(&n, matrix.data(), &inc)) + matrix.exponent();
    };
    
    template<typename Value>
    void density_matrix(Matrix<Host, Value>& dest, Matrix<Host, Value> const& B, Vector<Host> const& prop, Matrix<Host, Value> const& A, Energies<Host> const& energies, itf::Batcher<Value>& batcher) {
        dest.I() = A.I(); dest.J() = B.I(); dest.exponent() = A.exponent() + B.exponent();
        char conjNo = 'n'; char conjYes = 'c'; Value one = 1.; Value zero = .0;
        gemm(&conjYes, &conjNo, &B.I(), &A.I(), &A.J(), &one, B.data(), &B.J(), A.data(), &A.J(), &zero, dest.data(), &dest.J());
        
        dest.exponent() += prop.exponent(); double const deltaTime = -prop.time(); // this is confusing, change time -> -time
        auto d = dest.data(); auto const p = prop.data(); auto const e = energies.data();
        for(int i = 0; i < dest.I(); ++i)
            for(int j = 0; j < dest.J(); ++j) {
                double const deltaE = e[j] - e[i]; double const delta = deltaTime*deltaE;
                d[j + dest.J()*i] *= std::abs(delta) > 1.e-7 ? (p[i] - p[j])/deltaE : deltaTime/2.*(p[i] + p[j] + delta/2.*(p[j] - p[i])); //approximation is symmetric
            }
    };
    
    template<typename Value>
    void add(Matrix<Host, Value>& dest, ut::Zahl<Value> const& fact, Matrix<Host, Value> const& source, itf::Batcher<Value>& batcher) {
        int const one = 1; int const n = source.I()*source.J();
        Value const x = (fact*ut::Zahl<Value>(1., source.exponent())).get();
        axpy(&n, &x, source.data(), &one, dest.data(), &one);
    };
    
    template<typename Value>
    void add(Value* dest, Value fact, Matrix<Host, Value> const& source) {
        int const one = 1; int const n = source.I()*source.J();
        axpy(&n, &fact, source.data(), &one, dest, &one);
    };
}


#endif

