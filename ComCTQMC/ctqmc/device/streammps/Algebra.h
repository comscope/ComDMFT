#ifndef CTQMC_DEVICE_STREAMMPS_ALGEBRA_H
#define CTQMC_DEVICE_STREAMMPS_ALGEBRA_H

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

#include "../include/Allocator.h"

#include "../../include/Utilities.h"
#include "../../include/impurity/Algebra.h"

#include "../../../include/BlasLapack.h"
#include "../../../include/JsonX.h"
#include "../../../include/io/Matrix.h"


namespace imp {
    
    struct Device {};
    
    
    
    
    int pci_id_size();
    int get_pci_ids(std::vector<char>&);
    void init_device(std::vector<char> const&, std::size_t);
    void release_device();
    

    
    template<>
    struct Batcher<Device, double> : itf::Batcher<double> {
        enum class Phase { record, execute };
        
        Batcher() = delete;
        Batcher(std::size_t);
        Batcher(Batcher const&) = delete;
        Batcher(Batcher&&) = delete;
        Batcher& operator=(Batcher const&) = delete;
        Batcher& operator=(Batcher&&) = delete;
        
        double* get_callback(std::function<void(double)> callBack);
        
        int is_ready();
        void launch();
        
        ~Batcher();
    private:
        std::size_t size_;
        
        Phase phase_;
        
        std::vector<std::function<void(double)>>  callBack_;
        double* hostCallBackBuffer_;
        device::Allocator::Data<double> deviceCallBackBuffer_;
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
        Energies<Device> const& energies() const { return energies_;}; //????   pass eig to copyEvolveL etc ??
        
    private:
        double const time_;
        double const exponent_;
        Energies<Device> const& energies_;     //????   pass eig to copyEvolveL etc ??
    };
    
    
    template<>
    struct Matrix<Device, double> {
        using data_type = typename device::Allocator::template Data<double>;
        
        struct Identity { Identity(int d) : dim(d) {}; int const dim;};
        struct Zero { Zero(int d) : dim(d) {}; int const dim;};
        
        Matrix() = delete;
        Matrix(int size);
        Matrix(Identity const& identity);
        Matrix(Zero const& zero);
        Matrix(int I, int J, io::rmat const& mat);
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

    
    
    void copyEvolveL(Matrix<Device, double>& dest, Vector<Device> const& prop, Matrix<Device, double> const& source, itf::Batcher<double>& batcher);
    void mult(Matrix<Device, double>& dest, Matrix<Device, double> const& L, Matrix<Device, double> const& R, itf::Batcher<double>& batcher);
    void evolveL(Vector<Device> const& prop, Matrix<Device, double>& arg, itf::Batcher<double>& batcher);

    void trace(ut::Zahl<double>* Z, ut::Zahl<double>* accZ, Matrix<Device, double> const& matrix, itf::Batcher<double>& batcher);
    void traceAtB(ut::Zahl<double>* Z, ut::Zahl<double>* accZ, Matrix<Device, double> const& A, Matrix<Device, double> const& B, itf::Batcher<double>& batcher);
    void norm(double* norm, Matrix<Device, double> const& matrix, itf::Batcher<double>& batcher);

    void density_matrix(Matrix<Device, double>& dest, Matrix<Device, double> const& B, Vector<Device> const& prop, Matrix<Device, double> const& A, Energies<Device> const& energies, itf::Batcher<double>& batcher);  // async
    void add(Matrix<Device, double>& dest, ut::Zahl<double> const& fact, Matrix<Device, double> const& source, itf::Batcher<double>& batcher);

    void add(double* dest, double fact, Matrix<Device, double> const& source);    
    
}


#endif

