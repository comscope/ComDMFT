#ifndef CTQMC_DEVICE_INCLUDE_ERRCHK_H
#define CTQMC_DEVICE_INCLUDE_ERRCHK_H

#include <cuda.h>
#include <cuda_runtime_api.h>


#define cudaErrchk(ans) { cudaAssert((ans), __FILE__, __LINE__); }
    inline void cudaAssert(cudaError_t code, const char *file, int line) {
        if(code != cudaSuccess) throw std::runtime_error("cuda error: " + std::string(cudaGetErrorString(code)) + " in " + std::string(file) + " line " + std::to_string(line));
    }


#ifdef HAVE_CUBLAS

#include <cublas_v2.h>

#define cublasErrchk(ans) { cublasAssert((ans), __FILE__, __LINE__); }
inline void cublasAssert(cublasStatus_t code, const char *file, int line) {
    if(code != CUBLAS_STATUS_SUCCESS) throw std::runtime_error("cublas error");
}

#endif  //HAVE_CUBLAS


#endif  //CTQMC_DEVICE_INCLUDE_ERRCHK_H

