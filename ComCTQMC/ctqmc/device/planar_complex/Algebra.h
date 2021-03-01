#ifndef CTQMC_DEVICE_DP_ALGEBRA_H
#define CTQMC_DEVICE_DP_ALGEBRA_H

#include <stdexcept>
#include <list>
#include <set>
#include <functional>
#include <cassert>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <limits>
#include <vector>
#include <algorithm>

#include <thrust/complex.h>
#include <cuComplex.h>

#include "../include/Allocator.h"

#include "../../include/Utilities.h"
#include "../../include/impurity/Algebra.h"

#include "../../../include/BlasLapack.h"
#include "../../../include/JsonX.h"
#include "../../../include/io/Matrix.h"

namespace imp {


    //need an interface between std::complex<double> and thrust::complex<double> & cuDoubleComplex
    template <typename Value>
    struct cuda_value{};

    template <>
    struct cuda_value<double>{
        using type = double;
        using scalar = double;
        using cscalar = double;
        static constexpr std::size_t size = 1;
    };

    template <>
    struct cuda_value<ut::complex>{
        using type = double;
        using scalar = thrust::complex<double>;
        using cscalar = cuDoubleComplex; //CUDA doesn't always play nice with C++...
        static constexpr std::size_t size = 2;
        //planar complex : array of real followed by array of imaginary (both of type type)
    };

    template <typename Value> using cuda_value_t = typename cuda_value<Value>::type; //How the matrix is stored
    template <typename Value> using cuda_value_scalar = typename cuda_value<Value>::scalar; //How a scalar is stored (C++type)
    template <typename Value> using cuda_value_cscalar = typename cuda_value<Value>::cscalar; //How a scalar is stored (Ctype)

    
    struct Device {};
    
    
    
    
    int pci_id_size();
    int get_pci_ids(std::vector<char>&);
    void init_device(int const, std::size_t);
    void release_device();
    

    template <typename Value>
    struct Kernel;
    
    template<typename Value>
    struct Batcher<Device, Value> : itf::Batcher<Value> {
        enum class Phase { record, execute, finalize };
        
        Batcher() = delete;
        Batcher(std::size_t);
        Batcher(Batcher const&) = delete;
        Batcher(Batcher&&) = delete;
        Batcher& operator=(Batcher const&) = delete;
        Batcher& operator=(Batcher&&) = delete;
        ~Batcher();
        
        //generic get function
        cuda_value_scalar<Value>* get_callback(std::function<void(cuda_value_scalar<Value>)> callBack);
        
        template<typename K> K& get_kernel();
        
        int is_ready();
        void launch();
    
    private:
        std::size_t size_;
        std::size_t numberOfKernels_;
 
        Phase phase_;
        cudaStream_t stream_;

        Kernel<Value>* hostKernelBuffer_;
        device::Allocator::Data<Kernel<Value>> deviceKernelBuffer_;
        
        std::vector<std::function<void(cuda_value_scalar<Value>)>>  callBack_;
        cuda_value_scalar<Value>* hostCallBackBuffer_;
        device::Allocator::Data<cuda_value_scalar<Value>> deviceCallBackBuffer_;
        
        device::Allocator::Data<device::Byte> memory_;
        
    };
    
    
    
    
    template<>
    struct Energies<Device> {
        using data_type = typename device::Allocator::template Data<double>; //An Haesslichkeit schwer z'uebertreffe ......
        
        Energies() = delete;
        Energies(jsx::value const& jParams, std::vector<double> const& eig);
        Energies(Energies const&) = delete;
        Energies(Energies&&) = delete;
        Energies& operator=(Energies const&) = delete;
        Energies& operator=(Energies&&) = delete;
        ~Energies();
        
        int const& dim0() const { return dim0_;};
        int const& dim() const { return dim_;}
        double const& ln_dim() const { return ln_dim_;};
        data_type& data() { return data_;}
        data_type const& data() const { return data_;}
        double const& min() const { return min_;}
        
    private:
        int const dim0_;
        int const dim_;
        double const ln_dim_;
        data_type data_;
        double min_;
    };
    
    
    template<>
    struct Vector<Device> {
        Vector() = delete;
        Vector(double time, Energies<Device> const& energies);
        Vector(Vector const&) = delete;
        Vector(Vector&&) = delete;
        Vector& operator=(Vector const&) = delete;
        Vector& operator=(Vector&&) = delete;
        ~Vector();
        
        double const& time() const { return time_;};
        double const& exponent() const { return exponent_;}
        Energies<Device> const& energies() const { return energies_;}; //????   pass eig to copyEvolveLL etc ??
        
    private:
        double const time_;
        double const exponent_;
        Energies<Device> const& energies_;     //????   pass eig to copyEvolveLL etc ??
    };
    
    
    template<typename Value>
    struct Matrix<Device, Value> {
        using data_type = typename device::Allocator::template Data<cuda_value_t<Value>>;
        
        struct Identity { Identity(int d) : dim(d) {}; int const dim;};
        struct Zero { Zero(int d) : dim(d) {}; int const dim;};
        
        Matrix() = delete;
        Matrix(int size);
        Matrix(Identity const& identity);
        Matrix(Zero const& zero);
        Matrix(int I, int J, io::Matrix<Value> const& mat);
        
        Matrix(Matrix const&) = delete;
        Matrix(Matrix&&) = delete;
        Matrix& operator=(Matrix const&) = delete;
        Matrix& operator=(Matrix&&) = delete;
        ~Matrix();
        
        int& I() { return I_;}
        int& J() { return J_;}
        int const& I() const { return I_;}
        int const& J() const { return J_;}
        data_type& data() { return data_;}
        data_type const& data() const { return data_;}
        double& exponent() { return exponent_;}
        double const& exponent() const { return exponent_;}
        
    protected:        
        int I_, J_;
        data_type data_;
        double exponent_;
    };

    template <typename Value>
    void copyEvolveL(Matrix<Device, Value>& dest, Vector<Device> const& prop, Matrix<Device, Value> const& source, itf::Batcher<Value>& batcher);

    template <typename Value>
    void mult(Matrix<Device, Value>& dest, Matrix<Device, Value> const& L, Matrix<Device, Value> const& R, itf::Batcher<Value>& batcher);

    template <typename Value>
    void evolveL(Vector<Device> const& prop, Matrix<Device, Value>& arg, itf::Batcher<Value>& batcher);

    template <typename Value>
    void trace(ut::Zahl<Value>* Z, ut::Zahl<Value>* accZ, Matrix<Device, Value> const& matrix, itf::Batcher<Value>& batcher);

    template <typename Value>
    void traceAtB(ut::Zahl<Value>* Z, ut::Zahl<Value>* accZ, Matrix<Device, Value> const& A, Matrix<Device, Value> const& B, itf::Batcher<Value>& batcher);

    template <typename Value>
    void norm(double* norm, Matrix<Device, Value> const& matrix, itf::Batcher<Value>& batcher);

    template <typename Value>
    void density_matrix(Matrix<Device, Value>& dest, Matrix<Device, Value> const& B, Vector<Device> const& prop, Matrix<Device, Value> const& A, Energies<Device> const& energies, itf::Batcher<Value>& batcher);  // async

    template <typename Value>
    void add(Matrix<Device, Value>& dest, ut::Zahl<Value> const& fact, Matrix<Device, Value> const& source, itf::Batcher<Value>& batcher);

    template <typename Value>
    void add(Value* dest, Value fact, Matrix<Device, Value> const& source);



}


#endif

