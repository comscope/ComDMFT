#include <cstring>
#include <type_traits>
#include <sstream>
#include <stdexcept>

#include "Algebra.h"

#include "../include/Errchk.h"

#include "../../../include/BlasLapack.h"
#include "../../../include/mpi/Basic.h"


#ifdef HAVE_CUBLAS

cublasHandle_t handle;

#else

#include <cutlass/gemm/gemm.h>
#include <cutlass/gemm/dgemm_traits.h>

#endif


using namespace imp;
using namespace device;

constexpr int WarpSize = 32;
Allocator* alloc = nullptr;


int imp::pci_id_size() {
    return 16;
}

// @#$!#$@!# cudaGetDevicePciBusId is not working properly on SummitDev ......
void get_pci_id(char* pci_id, int deviceId) {
    std::stringstream temp;
    cudaDeviceProp deviceProperties; cudaErrchk(cudaGetDeviceProperties(&deviceProperties, deviceId));
    temp << std::hex << deviceProperties.pciDomainID << ":" << deviceProperties.pciBusID << ":" << deviceProperties.pciDeviceID;
    std::string str = temp.str(); std::copy(str.begin(), str.end(), pci_id);
}


int imp::get_pci_ids(std::vector<char>& pciIds) {
    int deviceCount; cudaErrchk(cudaGetDeviceCount(&deviceCount));
    
    pciIds.resize(deviceCount*pci_id_size(), '\0');
    for(int id = 0; id < deviceCount; ++id) get_pci_id(&pciIds[id*pci_id_size()], id);
    
    return deviceCount;
}

void imp::init_device(std::vector<char> const& pciId, std::size_t processesPerDevice) {
    int deviceId; cudaErrchk(cudaDeviceGetByPCIBusId(&deviceId, pciId.data()));
    cudaErrchk(cudaSetDevice(deviceId));
    
    cudaDeviceProp deviceProperties; cudaErrchk(cudaGetDeviceProperties(&deviceProperties, deviceId));
    if(deviceProperties.computeMode != cudaComputeModeExclusive && deviceProperties.computeMode != cudaComputeModeExclusiveProcess)
        throw std::runtime_error("Please set GPU compute mode to \"cudaComputeModeExclusive\" or \"cudaComputeModeExclusiveProcess\"");
    if(deviceProperties.warpSize != WarpSize)
        throw std::runtime_error("Please set WarpSize in AlgebraDevice.cu to " + std::to_string(deviceProperties.warpSize));
    
    cudaErrchk(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));

#ifdef HAVE_CUBLAS
    
    cublasErrchk(cublasCreate(&handle));
    
#endif
    
    alloc = new Allocator((0.8/processesPerDevice)*deviceProperties.totalGlobalMem);
}

void imp::release_device() {
    delete alloc; alloc = nullptr;
    
#ifdef HAVE_CUBLAS
    
    cublasErrchk(cublasDestroy(handle));
    
#endif

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//-------------------------------------------------------------------------------------------------------------------------------------------------

imp::Energies<Device>::Energies(jsx::value const& jParams, std::vector<double> const& energies) :
dim0_(energies.size()),
dim_(jParams.is("trunc dim") ? std::min<int>(dim0_, jParams("trunc dim").int64()) : dim0_),
ln_dim_(std::log(dim_)),
data_(alloc->get<double>(dim_)),
min_(std::numeric_limits<double>::max()) {
    for(int i = 0; i < dim_; ++i)
        min_ = std::min(min_, energies[i]);
    
    cudaErrchk(cudaMemcpy(data_.ptr(), energies.data(), dim_*sizeof(double), cudaMemcpyHostToDevice));
}
imp::Energies<Device>::~Energies() {
    alloc->free(data_);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

imp::Vector<Device>::Vector(double time, Energies<Device> const& energies) :
time_(time),
exponent_(time*energies.min()),
energies_(energies) {
}
imp::Vector<Device>::~Vector() {
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

imp::Matrix<Device, double>::Matrix(int size) : data_(alloc->get<double>(size)) {
}
imp::Matrix<Device, double>::Matrix(Matrix<Device, double>::Identity const& identity) : I_(identity.dim), J_(identity.dim), data_(alloc->get<double>(I_*J_)), exponent_(.0) {
    double* temp = new double[I_*J_]; std::memset(temp, 0, I_*J_*sizeof(double));
    for(int i = 0; i < identity.dim; ++i) temp[i*(identity.dim + 1)] = 1.; //kack memset isch das allgemein für double's ?
    
    cudaErrchk(cudaMemcpy(data_.ptr(), temp, I_*J_*sizeof(double), cudaMemcpyHostToDevice));
    
    delete[] temp;
}
imp::Matrix<Device, double>::Matrix(Matrix<Device, double>::Zero const& zero) : I_(zero.dim), J_(zero.dim), data_(alloc->get<double>(I_*J_)), exponent_(.0) {
    double* temp = new double[I_*J_]; std::memset(temp, 0, I_*J_*sizeof(double));
    
    cudaErrchk(cudaMemcpy(data_.ptr(), temp, I_*J_*sizeof(double), cudaMemcpyHostToDevice));
    
    delete[] temp;
}
imp::Matrix<Device, double>::Matrix(int I, int J, io::rmat const& mat) : I_(I), J_(J), data_(alloc->get<double>(I_*J_)), exponent_(.0) {
    double* temp = new double[I_*J_];
    
    for(int i = 0; i < I; ++i)
        for(int j = 0; j < J; ++j)
            temp[j + J*i] = mat(i, j);
    
    cudaErrchk(cudaMemcpy(data_.ptr(), temp, I_*J_*sizeof(double), cudaMemcpyHostToDevice));
    
    delete[] temp;
}
imp::Matrix<Device, double>::~Matrix() {
    alloc->free(data_);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

void imp::add(double* dest, double fact, Matrix<Device, double> const& source) {
    int const N = source.I()*source.J(); int const one = 1;
    double* temp = new double[N];                              //Ja scheisse das isch beschisse, passiert aber nit oft.
    
    cudaErrchk(cudaMemcpy(temp, source.data().ptr(), N*sizeof(double), cudaMemcpyDeviceToHost));
    daxpy_(&N, &fact, temp, &one, dest, &one);
    
    delete[] temp;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//-------------------------------------------------------------------------------------------------------------------------------------------------

struct CopyEvolveL {
    double time;
    double shift;
    double const* energies;
    double const* source;
    double* dest;
};

__global__ void kerCopyEvolveL(CopyEvolveL args) {
    int const i = blockIdx.x; int const j = threadIdx.x;

    args.dest[j + blockDim.x*i] = exp(args.time*args.energies[i] - args.shift)*args.source[j + blockDim.x*i];
};

void imp::copyEvolveL(Matrix<Device, double>& dest, Vector<Device> const& prop, Matrix<Device, double> const& source, itf::Batcher<double>& batcher) {
    dest.I() = source.I(); dest.J() = source.J(); dest.exponent() = source.exponent() + prop.exponent(); // eigentli source.exponent() = 0, isch aber sicherer so
    
    CopyEvolveL args;
    
    args.time     = prop.time();
    args.shift    = prop.exponent();
    args.energies = prop.energies().data().ptr();
    args.source   = source.data().ptr();
    args.dest     = dest.data().ptr();
    
    kerCopyEvolveL<<<source.I(), source.J()>>>(args);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

#ifndef HAVE_CUBLAS

void cutlassDgemm(
                         int M,
                         int N,
                         int K,
                         double const *A,
                         int lda,
                         double const *B,
                         int ldb,
                         double *C,
                         int ldc) {
    
    typedef cutlass::gemm::DgemmTraits<
    cutlass::MatrixLayout::kColumnMajor,
    cutlass::MatrixLayout::kColumnMajor
    >
    GemmTraits;
    
    typedef cutlass::gemm::Gemm<GemmTraits> Gemm;
    
    typename Gemm::Params params;
    
    params.initialize(M, N, K, 1., A, lda, B, ldb, .0, C, ldc, C, ldc);
    
    Gemm::launch(params);
}

#endif

void imp::mult(Matrix<Device, double>& dest, Matrix<Device, double> const& L, Matrix<Device, double> const& R, itf::Batcher<double>& batcher) {
    dest.I() = L.I(); dest.J() = R.J(); dest.exponent() = L.exponent() + R.exponent();
    
#ifdef HAVE_CUBLAS
    
    double one = 1.; double zero = .0;
    
    cublasErrchk(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, R.J(), L.I(), L.J(), &one, R.data().ptr(), R.J(), L.data().ptr(), L.J(), &zero, dest.data().ptr(), R.J()));
    
#else
    
    cutlassDgemm(R.J(), L.I(), L.J(), R.data().ptr(), R.J(), L.data().ptr(), L.J(), dest.data().ptr(), R.J());
    
#endif
    
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

struct EvolveL {
    double time;
    double shift;
    double const* energies;
    double* arg;
};

__global__ void kerEvolveL(EvolveL args) {
    int const i = blockIdx.x; int const j = threadIdx.x;
   
    args.arg[j + blockDim.x*i] *= exp(args.time*args.energies[i] - args.shift);
};

void imp::evolveL(Vector<Device> const& prop, Matrix<Device, double>& arg, itf::Batcher<double>& batcher) {
    arg.exponent() += prop.exponent();
    
    EvolveL args;
    
    args.time     = prop.time();
    args.shift    = prop.exponent();
    args.energies = prop.energies().data().ptr();
    args.arg      = arg.data().ptr();
    
    kerEvolveL<<<arg.I(), arg.J()>>>(args);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

#if __CUDACC_VER_MAJOR__ >= 9

__device__ __forceinline__ void reduceWarp(int const tid, double* data, double* result) {
    double temp;
    temp = data[tid + 16]; __syncwarp();
    data[tid] += temp;     __syncwarp();
    temp = data[tid + 8];  __syncwarp();
    data[tid] += temp;     __syncwarp();
    temp = data[tid + 4];  __syncwarp();
    data[tid] += temp;     __syncwarp();
    temp = data[tid + 2];  __syncwarp();
    data[tid] += temp;     __syncwarp();
    temp = data[tid + 1];  __syncwarp();
    data[tid] += temp;     __syncwarp();
    if(tid == 0) *result = *data;
};

#else

__device__ __forceinline__ void reduceWarp(int const tid, double volatile* data, double* result) {
    data[tid] += data[tid + 16];
    data[tid] += data[tid + 8];
    data[tid] += data[tid + 4];
    data[tid] += data[tid + 2];
    data[tid] += data[tid + 1];
    if(tid == 0) *result = *data;
};

#endif

template<int Size>
__device__ __forceinline__ void reduce(int const tid, double* data, double* result) {
    if(tid < Size/2) data[tid] += data[tid + Size/2];
    __syncthreads();
    reduce<Size/2>(tid, data, result);
};

template<>
__device__ __forceinline__ void reduce<WarpSize>(int const tid, double* data, double* result) {
    if(tid < WarpSize) reduceWarp(tid, data, result);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------

struct Trace {
    double const* arg;
    double* result;
    int dim;
};

template<int BlockDim>
__global__ void kerTrace(Trace args) {
    __shared__ double cache[BlockDim + 16];   // I do not want some threads in the reduceWarp to read stuff outside the cache ... nobody of the nVidia freaks seems to care about this (and probabely they are right) but I do not see why
    cache[threadIdx.x] = .0;
    int i = threadIdx.x;
    
    while(i < args.dim) {
        cache[threadIdx.x] += args.arg[(args.dim + 1)*i];
        
        i += BlockDim;
    }
    __syncthreads();
    
    reduce<BlockDim>(threadIdx.x, cache, args.result);
};

void imp::trace(ut::Zahl<double>* Z, ut::Zahl<double>* accZ, Matrix<Device, double> const& matrix, itf::Batcher<double>& batcher) {
    double exponent = matrix.exponent();
    
    Trace args;
    
    args.arg    = matrix.data().ptr();
    args.result = get<Device>(batcher).get_callback([=](double buffer) { ut::Zahl<double> temp(buffer, exponent); if(Z) *Z = temp; if(accZ) *accZ += temp;});
    args.dim    = matrix.I();
    
    kerTrace<WarpSize><<<1, WarpSize>>>(args);
}


//-------------------------------------------------------------------------------------------------------------------------------------------------

struct TraceAtB {
    double const* At;
    double const* B;
    double* result;
    int size;
};

template<int BlockDim>
__global__ void kerTraceAtB(TraceAtB args) {
    __shared__ double cache[BlockDim + 16];
    cache[threadIdx.x] = .0;
    int i = threadIdx.x;
    
    while(i < args.size) {
        cache[threadIdx.x] += args.At[i]*args.B[i];
        
        i += BlockDim;
    }
    __syncthreads();
    
    reduce<BlockDim>(threadIdx.x, cache, args.result);
};

void imp::traceAtB(ut::Zahl<double>* Z, ut::Zahl<double>* accZ, Matrix<Device, double> const& At, Matrix<Device, double> const& B, itf::Batcher<double>& batcher) {
    double exponent = At.exponent() + B.exponent();
    
    TraceAtB args;
    
    args.At     = At.data().ptr();
    args.B      = B.data().ptr();
    args.result = get<Device>(batcher).get_callback([=](double buffer) { ut::Zahl<double> temp(buffer, exponent); if(Z) *Z = temp; if(accZ) *accZ += temp;});
    args.size   = At.I()*At.J();
    
    kerTraceAtB<1024><<<1, 1024>>>(args);
}


//-------------------------------------------------------------------------------------------------------------------------------------------------

struct Norm {
    double const* arg;
    double* result;
    int size;
};

template<int BlockDim>
__global__ void kerNorm(Norm args) {
    __shared__ double cache[BlockDim + 16];
    cache[threadIdx.x] = .0;
    double value; int i = threadIdx.x;
    
    while(i < args.size) {
        value = args.arg[i];
        cache[threadIdx.x] += value*value;
        
        i += BlockDim;
    }
    __syncthreads();
    
    reduce<BlockDim>(threadIdx.x, cache, args.result);
};

void imp::norm(double* norm, Matrix<Device, double> const& matrix, itf::Batcher<double>& batcher) {
    double exponent = matrix.exponent();
    
    Norm args;
    
    args.arg    = matrix.data().ptr();
    args.result = get<Device>(batcher).get_callback([=](double buffer) { *norm = std::log(buffer)/2. + exponent;});
    args.size   = matrix.I()*matrix.J();
    
    kerNorm<1024><<<1, 1024>>>(args);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

void imp::density_matrix(Matrix<Device, double>& dest, Matrix<Device, double> const& B, Vector<Device> const& prop, Matrix<Device, double> const& A, Energies<Device> const& energies, itf::Batcher<double>& batcher) 
{
    throw std::runtime_error("imp::density_matrix: not implemented !");
};

//-------------------------------------------------------------------------------------------------------------------------------------------------

#ifndef HAVE_CUBLAS

__global__ void kerAdd(int size, double fact, double const* source, double* dest)
{
    int const index = blockDim.x*blockIdx.x + threadIdx.x;
    
    if(index < size) dest[index] += fact*source[index];
};

#endif


void imp::add(Matrix<Device, double>& dest, ut::Zahl<double> const& fact, Matrix<Device, double> const& source, itf::Batcher<double>& batcher)
{
    double alpha = (fact*ut::exp(source.exponent())).get();
    
#ifdef HAVE_CUBLAS
    
    cublasErrchk(cublasDaxpy(handle, source.I()*source.J(), &alpha, source.data().ptr(), 1, dest.data().ptr(), 1));
    
#else
    
    kerAdd<<<(source.I()*source.J() + 256 - 1)/256, 256>>>(source.I()*source.J(), alpha, source.data().ptr(), dest.data().ptr());
    
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

imp::Batcher<Device, double>::Batcher(std::size_t size) :
size_(size),
phase_(Phase::record),
deviceCallBackBuffer_(alloc->get<double>(size_)) {
    cudaErrchk(cudaMallocHost(reinterpret_cast<void**>(&hostCallBackBuffer_), size_*sizeof(double)));
}

double* imp::Batcher<Device, double>::get_callback(std::function<void(double)> callBack) {
    int index = callBack_.size();  callBack_.push_back(callBack);
    return deviceCallBackBuffer_.ptr() + index;
};

void imp::Batcher<Device, double>::launch() {
    if(callBack_.size()) {
        cudaErrchk(cudaMemcpy(hostCallBackBuffer_, deviceCallBackBuffer_.ptr(), callBack_.size()*sizeof(double), cudaMemcpyDeviceToHost));
        for(std::size_t index = 0; index < callBack_.size(); ++index) callBack_[index](hostCallBackBuffer_[index]);
        callBack_.clear();
    }
};

int imp::Batcher<Device, double>::is_ready() {
    if(callBack_.size()) throw std::runtime_error("Batcher::is_ready");
    
    return 1;
};

imp::Batcher<Device, double>::~Batcher() {
    alloc->free(deviceCallBackBuffer_);
    cudaErrchk(cudaFreeHost(hostCallBackBuffer_));
};






