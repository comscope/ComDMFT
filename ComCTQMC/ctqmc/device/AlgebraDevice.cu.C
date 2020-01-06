#include <cstring>
#include <type_traits>
#include <sstream>

#include "AlgebraDevice.h"
#include "Variant.h"
#include "../../include/BlasLapack.h"
#include "../../include/mpi/Basic.h"

using namespace imp;

constexpr int WarpSize = 32;
int deviceId;
AllocDevice* alloc = nullptr;

cublasHandle_t hostHandle;

__device__ double deviceZero = .0;
__device__ double deviceOne  = 1.;
__device__ cublasHandle_t deviceHandle;

__global__ void kerCublasHandle(int const flag) {
    flag ? cublasCreate(&deviceHandle) : cublasDestroy(deviceHandle);
};

int imp::pci_id_size() {
    return 16;
}

// @#$!#$@!# cudaGetDevicePciBusId is not working properly on SummitDev ......
void get_pci_id(char* pci_id, int deviceId) {
    std::stringstream stream;
    cudaDeviceProp deviceProperties; cudaGetDeviceProperties(&deviceProperties, deviceId);
    stream << std::hex << deviceProperties.pciDomainID << ":" << deviceProperties.pciBusID << ":" << deviceProperties.pciDeviceID;
    std::string str = stream.str(); std::copy(str.begin(), str.end(), pci_id);
}


int imp::get_pci_ids(std::vector<char>& pciIds) {
    int deviceCount; cudaGetDeviceCount(&deviceCount);
    
    pciIds.resize(deviceCount*pci_id_size(), '\0');
    for(int id = 0; id < deviceCount; ++id) get_pci_id(&pciIds[id*pci_id_size()], id);

    return deviceCount;
}

void imp::init_device(std::vector<char> const& pciId) {
    cudaDeviceGetByPCIBusId(&deviceId, pciId.data()); cudaSetDevice(deviceId);
    
    cudaDeviceProp deviceProperties; cudaGetDeviceProperties(&deviceProperties, deviceId);
    if(deviceProperties.computeMode != cudaComputeModeExclusive && deviceProperties.computeMode != cudaComputeModeExclusiveProcess)
        throw std::runtime_error("Please set GPU compute mode to \"cudaComputeModeExclusive\" or \"cudaComputeModeExclusiveProcess\"");
    if(deviceProperties.warpSize != WarpSize)
        throw std::runtime_error("Please set WarpSize in AlgebraDevice.cu to " + std::to_string(deviceProperties.warpSize));

    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
            
    cublasCreate(&hostHandle);
    kerCublasHandle<<<1, 1>>>(1);
    cudaDeviceSynchronize();
    
    alloc = new AllocDevice(0.8*deviceProperties.totalGlobalMem);
}

void imp::release_device() {
    delete alloc; alloc = nullptr;
    kerCublasHandle<<<1, 1>>>(0);
    cublasDestroy(hostHandle);
    cudaDeviceSynchronize();
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<>
imp::Energies<AllocDevice>::Energies(jsx::value const& jParams, std::vector<double> const& energies) :
dim0_(energies.size()),
dim_(jParams.is("trunc dim") ? std::min<int>(dim0_, jParams("trunc dim").int64()) : dim0_),
ln_dim_(std::log(dim_)),
data_(alloc->get<double>(dim_)),
min_(std::numeric_limits<double>::max()) {
    for(int i = 0; i < dim_; ++i)
        min_ = std::min(min_, energies[i]);

    cudaMemcpy(data_.ptr(), energies.data(), dim_*sizeof(double), cudaMemcpyHostToDevice);
}
template<>
imp::Energies<AllocDevice>::~Energies() {
    alloc->free(data_);
}

template<>
imp::Vector<AllocDevice>::Vector(double time, Energies<AllocDevice> const& energies) :
time_(time),
exponent_(time*energies.min()),
energies_(energies),
data_(alloc->null_ptr<double>()) {
}
template<>
imp::Vector<AllocDevice>::~Vector() {
    alloc->free(data_);
}

template<>
imp::Matrix<AllocDevice>::Matrix(int size) : data_(alloc->get<double>(size)) {
}
template<>
imp::Matrix<AllocDevice>::Matrix(Matrix<AllocDevice>::Identity const& identity) : I_(identity.dim), J_(identity.dim), data_(alloc->get<double>(I_*J_)), exponent_(.0) {
    double* temp = new double[I_*J_]; std::memset(temp, 0, I_*J_*sizeof(double));
    for(int i = 0; i < identity.dim; ++i) temp[i*(identity.dim + 1)] = 1.; //kack memset isch das allgemein für double's ?
    
    cudaMemcpy(data_.ptr(), temp, I_*J_*sizeof(double), cudaMemcpyHostToDevice);
    
    delete[] temp;
}
template<>
imp::Matrix<AllocDevice>::Matrix(Matrix<AllocDevice>::Zero const& zero) : I_(zero.dim), J_(zero.dim), data_(alloc->get<double>(I_*J_)), exponent_(.0) {
    double* temp = new double[I_*J_]; std::memset(temp, 0, I_*J_*sizeof(double));
    
    cudaMemcpy(data_.ptr(), temp, I_*J_*sizeof(double), cudaMemcpyHostToDevice);
    
    delete[] temp;
}
template<>
imp::Matrix<AllocDevice>::Matrix(int I, int J, io::rmat const& mat) : I_(I), J_(J), data_(alloc->get<double>(I_*J_)), exponent_(.0) {
    double* temp = new double[I_*J_];
    
    for(int i = 0; i < I; ++i)
        for(int j = 0; j < J; ++j)
            temp[j + J*i] = mat(i, j);
    
    cudaMemcpy(data_.ptr(), temp, I_*J_*sizeof(double), cudaMemcpyHostToDevice);
    
    delete[] temp;
}
template<>
imp::Matrix<AllocDevice>::~Matrix() {
    alloc->free(data_);
}

void imp::axpy(double* dest, double fact, Matrix<AllocDevice> const& source) {
    int const N = source.I()*source.J(); int const one = 1;
    double* temp = new double[N];                              //Ja scheisse das isch beschisse, passiert aber nit oft.
    
    cudaMemcpy(temp, source.data().ptr(), N*sizeof(double), cudaMemcpyDeviceToHost);
    daxpy_(&N, &fact, temp, &one, dest, &one);
    
    delete[] temp;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void kerCopyEvolveL(double time, double const* energies, double shift, double* dest, double const* source) {
    int const i = blockIdx.x; int const j = threadIdx.x;
    
    dest[j + blockDim.x*i] = exp(time*energies[i] - shift)*source[j + blockDim.x*i];
};

struct CopyEvolveL {
    double time;
    double shift;
    double const* energies;
    double const* source;
    double* dest;
    int I;
    int J;
};

void imp::copyEvolveL(Matrix<AllocDevice>& dest, Vector<AllocDevice> const& prop, Matrix<AllocDevice> const& source, Batcher<AllocDevice>& batcher) {
    dest.I() = source.I(); dest.J() = source.J(); dest.exponent() = source.exponent() + prop.exponent(); // eigentli source.exponent() = 0, isch aber sicherer so
    
    auto& args = batcher.get_kernel<CopyEvolveL>();
    
    args.time     = prop.time();
    args.shift    = prop.exponent();
    args.energies = prop.energies().data().ptr();
    args.source   = source.data().ptr();
    args.dest     = dest.data().ptr();
    args.I        = source.I();
    args.J        = source.J();
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

struct Mult {
    double const* A;
    double const* B;
    double* C;
    int M;
    int N;
    int K;
};

void imp::mult(Matrix<AllocDevice>& dest, Matrix<AllocDevice> const& L, Matrix<AllocDevice> const& R, Batcher<AllocDevice>& batcher) {
    dest.I() = L.I(); dest.J() = R.J(); dest.exponent() = L.exponent() + R.exponent();
    
    auto& args = batcher.get_kernel<Mult>();
    
    args.A = L.data().ptr();
    args.B = R.data().ptr();
    args.C = dest.data().ptr();
    args.M = L.I();
    args.N = R.J();
    args.K = L.J();
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

__global__ void kerEvolveL(double time, double const* energies, double shift, double* arg) {
    int const i = blockIdx.x; int const j = threadIdx.x;
    
    arg[j + blockDim.x*i] *= exp(time*energies[i] - shift);
};

struct EvolveL {
    double time;
    double shift;
    double const* energies;
    double* arg;
    int I;
    int J;
};

void imp::evolveL(Vector<AllocDevice> const& prop, Matrix<AllocDevice>& arg, Batcher<AllocDevice>& batcher) {
    arg.exponent() += prop.exponent();
    
    auto& args = batcher.get_kernel<EvolveL>();
    
    args.time     = prop.time();
    args.shift    = prop.exponent();
    args.energies = prop.energies().data().ptr();
    args.arg      = arg.data().ptr();
    args.I        = arg.I();
    args.J        = arg.J();
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


template<int BlockDim>
__global__ void kerTrace(int dim, double const* arg, double* result) {
    __shared__ double cache[BlockDim + 16];   // I do not want some threads in the reduceWarp to read stuff outside the cache ... nobody of the nVidia freaks seems to care about this (and probabely they are right) but I do not see why
    cache[threadIdx.x] = .0;
    int i = threadIdx.x;
    
    while(i < dim) {
        cache[threadIdx.x] += arg[(dim + 1)*i];
        
        i += BlockDim;
    }
    __syncthreads();
    
    reduce<BlockDim>(threadIdx.x, cache, result);
};


struct Trace {
    double const* arg;
    double* result;
    int dim;
};

void imp::trace(ut::Zahl* Z, ut::Zahl* accZ, Matrix<AllocDevice> const& matrix, Batcher<AllocDevice>& batcher) {
    auto& args = batcher.get_kernel<Trace>(); double exponent = matrix.exponent();
    
    args.arg    = matrix.data().ptr();
    args.result = batcher.get_callback([=](double buffer) { ut::Zahl temp(buffer, exponent); if(Z) *Z = temp; if(accZ) *accZ += temp;});
    args.dim    = matrix.I();
}


//-------------------------------------------------------------------------------------------------------------------------------------------------

template<int BlockDim>
__global__ void kerTraceAtB(int const size, double const* At, double const* B, double* result) {
    __shared__ double cache[BlockDim + 16];
    cache[threadIdx.x] = .0;
    double valueAt, valueB; int i = threadIdx.x;
    
    while(i < size) {
        valueAt = At[i]; valueB = B[i];
        cache[threadIdx.x] += valueAt*valueB;
        
        i += BlockDim;
    }
    __syncthreads();
    
    reduce<BlockDim>(threadIdx.x, cache, result);
};

struct TraceAtB {
    double const* At;
    double const* B;
    double* result;
    int size;
};

void imp::traceAtB(ut::Zahl* Z, ut::Zahl* accZ, Matrix<AllocDevice> const& At, Matrix<AllocDevice> const& B, Batcher<AllocDevice>& batcher) {
    auto& args = batcher.get_kernel<TraceAtB>(); double exponent = At.exponent() + B.exponent();
    
    args.At     = At.data().ptr();
    args.B      = B.data().ptr();
    args.result = batcher.get_callback([=](double buffer) { ut::Zahl temp(buffer, exponent); if(Z) *Z = temp; if(accZ) *accZ += temp;});
    args.size   = At.I()*At.J();
}


//-------------------------------------------------------------------------------------------------------------------------------------------------

template<int BlockDim>
__global__ void kerNorm(int const size, double const* arg, double* result) {
    __shared__ double cache[BlockDim + 16];
    cache[threadIdx.x] = .0;
    double value; int i = threadIdx.x;
    
    while(i < size) {
        value = arg[i];
        cache[threadIdx.x] += value*value;
        
        i += BlockDim;
    }
    __syncthreads();
    
    reduce<BlockDim>(threadIdx.x, cache, result);
};

struct Norm {
    double const* arg;
    double* result;
    int size;
};

void imp::norm(double* norm, Matrix<AllocDevice> const& matrix, Batcher<AllocDevice>& batcher) {
    auto& args = batcher.get_kernel<Norm>(); double exponent = matrix.exponent();
    
    args.arg    = matrix.data().ptr();
    args.result = batcher.get_callback([=](double buffer) { *norm = std::log(buffer)/2. + exponent;});
    args.size   = matrix.I()*matrix.J();
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

void imp::density_matrix(Matrix<AllocDevice>& dest, ut::Zahl const& fact, Matrix<AllocDevice> const& Bt, Vector<AllocDevice> const& prop, Matrix<AllocDevice> const& A, Matrix<AllocDevice>& buffer, Batcher<AllocDevice>& batcher) {
    throw std::runtime_error("imp::density_matrix: not implemented !");
};

//-------------------------------------------------------------------------------------------------------------------------------------------------

struct Axpy {
    double const* source;
    double* dest;
    double fact;
    int size;
};

void imp::axpy(Matrix<AllocDevice>& dest, ut::Zahl const& fact, Matrix<AllocDevice> const& source, Batcher<AllocDevice>& batcher) {
    auto& args = batcher.get_kernel<Axpy>();
    
    args.source    = source.data().ptr();
    args.dest      = dest.data().ptr();
    args.fact      = (fact*ut::exp(source.exponent())).to_double();
    args.size      = source.I()*source.J();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef var::variant<CopyEvolveL, Mult, EvolveL, Trace, TraceAtB, Norm, Axpy> KerArgs;

struct alignas(16) imp::Kernel {
    KerArgs args;
    int id;
};

__global__ void kerLauncher(Kernel* kernel, int const N) {
    for(int n = 0; n < N; ++n) {
        Kernel ker = kernel[n];
        
        if(ker.id == var::index<Mult, KerArgs>::value) {
            
            auto& args = var::device::get<Mult>(ker.args);
            cublasDgemm(deviceHandle, CUBLAS_OP_N, CUBLAS_OP_N, args.N, args.M, args.K, &deviceOne, args.B, args.N, args.A, args.K, &deviceZero, args.C, args.N);
            
        } else if(ker.id == var::index<Norm, KerArgs>::value) {
            
            auto& args = var::device::get<Norm>(ker.args);
            kerNorm<1024><<<1, 1024>>>(args.size, args.arg, args.result);
            
        } else if(ker.id == var::index<EvolveL, KerArgs>::value) {
            
            auto& args = var::device::get<EvolveL>(ker.args);
            kerEvolveL<<<args.I, args.J>>>(args.time, args.energies, args.shift, args.arg);
            
        } else if(ker.id == var::index<CopyEvolveL, KerArgs>::value) {
            
            auto& args = var::device::get<CopyEvolveL>(ker.args);
            kerCopyEvolveL<<<args.I, args.J>>>(args.time, args.energies, args.shift, args.dest, args.source);
            
        } else if(ker.id == var::index<Trace, KerArgs>::value) {
            
            auto& args = var::device::get<Trace>(ker.args);
            kerTrace<WarpSize><<<1, WarpSize>>>(args.dim, args.arg, args.result);
            
        } else if(ker.id == var::index<TraceAtB, KerArgs>::value) {
            
            auto& args = var::device::get<TraceAtB>(ker.args);
            kerTraceAtB<1024><<<1, 1024>>>(args.size, args.At, args.B, args.result);
            
        } else if(ker.id == var::index<Axpy, KerArgs>::value) {
            
            auto& args = var::device::get<Axpy>(ker.args);
            cublasDaxpy(deviceHandle, args.size, &args.fact, args.source, 1, args.dest, 1);
            
        }
        
    }            
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

imp::Batcher<AllocDevice>::Batcher(std::size_t size) : 
size_(size), 
numberOfKernels_(0),
phase_(Phase::record),
deviceKernelBuffer_(alloc->get<Kernel>(size_)), 
deviceCallBackBuffer_(alloc->get<double>(size_)) {
    cudaStreamCreate(&stream_);

    cudaMallocHost(reinterpret_cast<void**>(&hostKernelBuffer_), size_*sizeof(Kernel));
    cudaMallocHost(reinterpret_cast<void**>(&hostCallBackBuffer_), size_*sizeof(double));
}

double* imp::Batcher<AllocDevice>::get_callback(std::function<void(double)> callBack) {
    if(phase_ != Phase::record) throw std::runtime_error("imp::Batcher::get_callback");
    
    int index = callBack_.size();  callBack_.push_back(callBack);
    return deviceCallBackBuffer_.ptr() + index;
};

template<typename K>
K& imp::Batcher<AllocDevice>::get_kernel() {
    if(phase_ != Phase::record) throw std::runtime_error("imp::Batcher::get_kernel");
    
    Kernel& ker = hostKernelBuffer_[numberOfKernels_++];
    ker.id = var::index<K, KerArgs>::value;
    return var::host::get<K>(ker.args);
};

void imp::Batcher<AllocDevice>::launch() {
    if(phase_ != Phase::record) throw std::runtime_error("imp::Batcher::launch");
    
    if(numberOfKernels_) {
        cudaMemcpyAsync(deviceKernelBuffer_.ptr(), hostKernelBuffer_, numberOfKernels_*sizeof(Kernel), cudaMemcpyHostToDevice, stream_);
        kerLauncher<<<1, 1, deviceId, stream_>>>(deviceKernelBuffer_.ptr(), numberOfKernels_);
        
        numberOfKernels_ = 0; phase_ = Phase::execute;
    }
};

int imp::Batcher<AllocDevice>::is_ready() {
    if(phase_ == Phase::execute) {
        if(cudaSuccess != cudaStreamQuery(stream_)) return 0;
        
        if(callBack_.size()) {
            cudaMemcpyAsync(hostCallBackBuffer_, deviceCallBackBuffer_.ptr(), callBack_.size()*sizeof(double), cudaMemcpyDeviceToHost, stream_);
            phase_ = Phase::finalize; return 0;
        }
    }

    if(phase_ == Phase::finalize) {
        if(cudaSuccess != cudaStreamQuery(stream_)) return 0;
        
        for(std::size_t index = 0; index < callBack_.size(); ++index) callBack_[index](hostCallBackBuffer_[index]);
        callBack_.clear(); 
    }

    phase_ = Phase::record; return 1;
};

imp::Batcher<AllocDevice>::~Batcher() {
    alloc->free(deviceCallBackBuffer_);
    cudaFreeHost(hostCallBackBuffer_);

    alloc->free(deviceKernelBuffer_);
    cudaFreeHost(hostKernelBuffer_);

    cudaStreamDestroy(stream_);
};





