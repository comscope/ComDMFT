#ifndef ALGEBRADEVICE
#define ALGEBRADEVICE

#include <stdexcept>
#include <list>
#include <set>
#include <cassert>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <limits>
#include <vector>
#include <algorithm>

#include "../include/impurity/Algebra.h"
#include "Allocator.h"

namespace imp {
    int pci_id_size();
    int get_pci_ids(std::vector<char>&);
    void init_device(std::vector<char> const&);
    void release_device();
    
    struct Kernel;
    
    template<>
    struct Batcher<AllocDevice> {
        enum class Phase { record, execute, finalize };
        
        Batcher() = delete;
        Batcher(std::size_t);
        Batcher(Batcher const&) = delete;
        Batcher(Batcher&&) = delete;
        Batcher& operator=(Batcher const&) = delete;
        Batcher& operator=(Batcher&&) = delete;
        
        double* get_callback(std::function<void(double)> callBack);
        template<typename K> K& get_kernel();
        
        int is_ready();
        void launch();
        
        ~Batcher();
    private:
        std::size_t size_;
        std::size_t numberOfKernels_;
 
        Phase phase_;
        cudaStream_t stream_;

        Kernel* hostKernelBuffer_;
        AllocDevice::Data<Kernel> deviceKernelBuffer_;
        
        std::vector<std::function<void(double)>>  callBack_;
        double* hostCallBackBuffer_;
        AllocDevice::Data<double> deviceCallBackBuffer_;
    };
    
    void copyEvolveL(Matrix<AllocDevice>& dest, Vector<AllocDevice> const& prop, Matrix<AllocDevice> const& source, Batcher<AllocDevice>& batcher);
    void mult(Matrix<AllocDevice>& dest, Matrix<AllocDevice> const& L, Matrix<AllocDevice> const& R, Batcher<AllocDevice>& batcher);
    void evolveL(Vector<AllocDevice> const& prop, Matrix<AllocDevice>& arg, Batcher<AllocDevice>& batcher);
    
    void trace(ut::Zahl* Z, ut::Zahl* accZ, Matrix<AllocDevice> const& matrix, Batcher<AllocDevice>& batcher);
    void traceAtB(ut::Zahl* Z, ut::Zahl* accZ, Matrix<AllocDevice> const& At, Matrix<AllocDevice> const& B, Batcher<AllocDevice>& batcher);
    void norm(double* norm, Matrix<AllocDevice> const& matrix, Batcher<AllocDevice>& batcher);
    
    /*not yet implemented*/ void density_matrix(Matrix<AllocDevice>& dest, ut::Zahl const& fact, Matrix<AllocDevice> const& Bt, Vector<AllocDevice> const& prop, Matrix<AllocDevice> const& A, Matrix<AllocDevice>& buffer, Batcher<AllocDevice>& batcher);  // async
    void axpy(Matrix<AllocDevice>& dest, ut::Zahl const& fact, Matrix<AllocDevice> const& source, Batcher<AllocDevice>& batcher);
    
    void axpy(double* dest, double fact, Matrix<AllocDevice> const& source);
}


#endif
