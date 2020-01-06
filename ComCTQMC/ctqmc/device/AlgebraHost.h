#ifndef ALGEBRAHOST_H
#define ALGEBRAHOST_H

#include <list>
#include <set>
#include <stdexcept>
#include <vector>

#include "../include/impurity/Algebra.h"

namespace imp {
    struct AllocHost {
        template<typename T> struct Data {
            Data() = delete;
            Data(std::size_t size) : ptr_(new T[size]) {};
            Data(Data const&) = delete;
            Data(Data&&) = delete;
            Data& operator=(Data const&) = delete;
            Data& operator=(Data&&) = delete;
            T* ptr() { return ptr_;}
            T const* ptr() const { return ptr_;}
            ~Data() { delete[] ptr_;};
        private:
            T* ptr_;
        };
    };
    
    template<>
    struct Batcher<AllocHost> {
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
    
    void copyEvolveL(Matrix<AllocHost>& dest, Vector<AllocHost> const& prop, Matrix<AllocHost> const& source, Batcher<AllocHost>& batcher);  //async
    void mult(Matrix<AllocHost>& dest, Matrix<AllocHost> const& L, Matrix<AllocHost> const& R, Batcher<AllocHost>& batcher);        //async
    void evolveL(Vector<AllocHost> const& prop, Matrix<AllocHost>& arg, Batcher<AllocHost>& batcher);                    //async
    
    void trace(ut::Zahl* Z, ut::Zahl* accZ, Matrix<AllocHost> const& matrix, Batcher<AllocHost>& batcher);    //async
    void traceAtB(ut::Zahl* Z, ut::Zahl* accZ, Matrix<AllocHost> const& At, Matrix<AllocHost> const& B, Batcher<AllocHost>& batcher);
    void norm(double* norm, Matrix<AllocHost> const& matrix, Batcher<AllocHost>& batcher);                       // async
    
    void density_matrix(Matrix<AllocHost>& dest, ut::Zahl const& fact, Matrix<AllocHost> const& Bt, Vector<AllocHost> const& prop, Matrix<AllocHost> const& A, Matrix<AllocHost>& buffer, Batcher<AllocHost>& batcher);  // async
    void axpy(Matrix<AllocHost>& dest, ut::Zahl const& fact, Matrix<AllocHost> const& source, Batcher<AllocHost>& batcher); // async
    
    void axpy(double* dest, double fact, Matrix<AllocHost> const& source);
}


#endif

