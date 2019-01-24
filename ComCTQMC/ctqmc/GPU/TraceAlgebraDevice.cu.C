#include <cstring>

#include "TraceAlgebraDevice.h"

#define SWITCH(arg) arg##Device
#include "../include/TraceAlgebra.h"
#undef SWITCH

using namespace TrDevice;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TrDevice::EigenValues::EigenValues(Json::Value const& jParams, Json::Value const& jBloc, double Ur0) :
dim0_(jBloc("Dimension").int64()),
dim_(jParams.is("TRUNC_DIM") ? std::min<int>(dim0_, jParams("TRUNC_DIM").int64()) : dim0_),
data_(Comm::alloc(dim_)),
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
    
    double* temp = new double[dim_]; // Nit noetig im prinzip .... 
    
    int const N = jBloc("Filling").int64();
    double const mu = jParams("mu").real64();
    for(int i = 0; i < dim_; ++i) {
        temp[i] = eig[i] - mu*N + .5*Ur0*N*N;  //eig[i] -(mu - .5*Ur0)*N + .5*Ur0*N*(N - 1)
        min_ = std::min(min_, data_[i]);
    }
    
    cudaMemcpy(get(data_), temp, dim_*sizeof(double), cudaMemcpyHostToDevice);
    
    delete[] temp;
}

TrDevice::EigenValues::~EigenValues() { Comm::free(data_);}

TrDevice::Vector::Vector(double time, EigenValues const& eig) :
time_(time),
eig_(eig),
exponent_(time*eig.min())  {
}

TrDevice::Vector::~Vector() {}

TrDevice::Matrix::Matrix(int size) : data_(Comm::alloc(size)) {
}

TrDevice::Matrix::Matrix(Matrix::Identity const& identity) : I_(identity.dim), J_(identity.dim), data_(Comm::alloc(I_*J_)), exponent_(.0) {
    double* temp = new double[I_*J_]; std::memset(temp, 0, I_*J_*sizeof(double));
    for(int i = 0; i < identity.dim; ++i) temp[i*(identity.dim + 1)] = 1.; //huere memset isch das allgemein für double's ?
    
    cudaMemcpy(get(data_), temp, I_*J_*sizeof(double), cudaMemcpyHostToDevice);
    
    delete[] temp;
}


TrDevice::Matrix::Matrix(Matrix::Zero const& zero) : I_(zero.dim), J_(zero.dim), data_(Comm::alloc(I_*J_)), exponent_(.0) {
    double* temp = new double[I_*J_]; std::memset(temp, 0, I_*J_*sizeof(double));
    
    cudaMemcpy(get(data_), temp, I_*J_*sizeof(double), cudaMemcpyHostToDevice);
    
    delete[] temp;
}

TrDevice::Matrix::Matrix(int I, int J, std::vector<double> const& mat, int dataColMajor, int I0, int J0) : I_(I), J_(J), data_(Comm::alloc(I_*J_)), exponent_(.0) {
    if(static_cast<int>(mat.size()) != I0*J0) throw(std::runtime_error("Wrong matrix size"));
    
    double* temp = new double[I_*J_];
    for(int i = 0; i < I; ++i) for(int j = 0; j < J; ++j) temp[j + J*i] = dataColMajor ? mat[i + j*I0] : mat[j + i*J0];
    
    cudaMemcpy(get(data_), temp, I_*J_*sizeof(double), cudaMemcpyHostToDevice);
    
    delete[] temp;
}

TrDevice::Matrix::~Matrix() {
    Comm::free(data_);
}

void TrDevice::axpy(double* dest, double fact, Matrix const& source) {
    int const N = source.I_*source.J_; int const one = 1;
    double* temp = new double[N];                              //Ja scheisse das isch beschisse, passiert aber nit oft.
    
    cudaMemcpy(temp, get(source.data_), N*sizeof(double), cudaMemcpyDeviceToHost);
    daxpy_(&N, &fact, temp, &one, dest, &one);
    
    delete[] temp;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int const WarpSize = 32;

__device__ double device_zero = .0;
__device__ double device_one  = 1.;
__device__ cublasHandle_t device_handle;


__global__ void kerExpH(int dim, double time, double const* eig, double shift, double* arg, int const inc = 1) {
    int const i = threadIdx.x + blockIdx.x*blockDim.x;
    
    if(i < dim) arg[inc*i] = exp(time*eig[i] - shift);
};


//-------------------------------------------------------------------------------------------------------------------------------------------------

struct MatrixExpH {
    static int const ID = 1;
    
    double time;
    double shift;
    double const* eig;
    double* arg;
    int dim;
};

//-------------------------------------------------------------------------------------------------------------------------------------------------

__global__ void kerCopyEvolveL(double time, double const* eig, double shift, double* dest, double const* source) {
    int const i = blockIdx.x; int const j = threadIdx.x;
    
    dest[j + blockDim.x*i] = exp(time*eig[i] - shift)*source[j + blockDim.x*i];
};

struct CopyEvolveL {
    static int const ID = 2;
    
    double time;
    double shift;
    double const* eig;
    double const* source;
    double* dest;
    int I;
    int J;
};

//-------------------------------------------------------------------------------------------------------------------------------------------------

struct Mult {
    static int const ID = 3;
    
    double const* A;
    double const* B;
    double* C;
    int M;
    int N;
    int K;
};

//-------------------------------------------------------------------------------------------------------------------------------------------------

__global__ void kerEvolveL(double time, double const* eig, double shift, double* arg) {
    int const i = blockIdx.x; int const j = threadIdx.x;
    
    arg[j + blockDim.x*i] *= exp(time*eig[i] - shift);
};

struct EvolveL {
    static int const ID = 4;
    
    double time;
    double shift;
    double const* eig;
    double* arg;
    int I;
    int J;
};

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<int _BlockDim>
__device__ void reduce(int const tid, double* data, double* result) {
    if(tid < _BlockDim/2) data[tid] = data[tid] + data[tid + _BlockDim/2];
    if(_BlockDim/2 > 32) { __syncthreads();};
    reduce<_BlockDim/2>(tid, data, result);
};

template<>
__device__ void reduce<1>(int const tid, double* data, double* result) {
    if(tid == 0) *result = *data;
};


template<int _BlockDim>
__global__ void kerTrace(int dim, double const* arg, double* result) {
    __shared__ double cache[_BlockDim];
    cache[threadIdx.x] = .0;
    int i = threadIdx.x;
    
    while(i < dim) {
        cache[threadIdx.x] += arg[(dim + 1)*i];
        
        i += _BlockDim;
    }
    __syncthreads();
    
    reduce<_BlockDim>(threadIdx.x, cache, result);
};

struct Trace {
    static int const ID = 5;
    
    double const* arg;
    double* result;
    int dim;
};

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<int _BlockDim>
__global__ void kerAccE(int dim, double const* eig, double const* arg, double* result) {
    __shared__ double cache[_BlockDim];
    cache[threadIdx.x] = .0;
    int i = threadIdx.x;
    
    while(i < dim) {
        cache[threadIdx.x] += eig[i]*arg[(dim + 1)*i];
        
        i += _BlockDim;
    }
    __syncthreads();
    
    reduce<_BlockDim>(threadIdx.x, cache, result);
};

struct AccE {
    static int const ID = 6;
    
    double const* eig;
    double const* arg;
    double* result;
    int dim;
};

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<int _BlockDim>
__global__ void kerNorm(int const size, double const* arg, double* result) {
    __shared__ double cache[_BlockDim];
    cache[threadIdx.x] = .0;
    double value; int i = threadIdx.x;
    
    while(i < size) {
        value = arg[i];
        cache[threadIdx.x] += value*value;
        
        i += _BlockDim;
    }
    __syncthreads();
    
    reduce<_BlockDim>(threadIdx.x, cache, result);
};

struct Norm {
    static int const ID = 7;
    
    double const* arg;
    double* result;
    int size;
};

struct Axpy {
    static int const ID = 8;
    
    double const* source;
    double* dest;
    double fact;
    int size;
};

//-------------------------------------------------------------------------------------------------------------------------------------------------

struct alignas(16) TrDevice::Kernel {
    union Get {
        MatrixExpH  matrixExpH;
        CopyEvolveL copyEvolveL;
        Mult        mult;
        EvolveL     evolveL;
        Trace       trace;
        AccE        accE;
        Norm        norm;
        Axpy        axpy;
    } get;
    
    int id;
};


__global__ void kerLauncher(Kernel* kernel, int const N) {
    for(int n = 0; n < N; ++n) {
        Kernel ker = kernel[n];
        
        if(ker.id == Mult::ID) {
            
            auto args = ker.get.mult;
            cublasDgemm(device_handle, CUBLAS_OP_N, CUBLAS_OP_N, args.N, args.M, args.K, &device_one, args.B, args.N, args.A, args.K, &device_zero, args.C, args.N);
            
        } else if(ker.id == Norm::ID) {
            
            auto args = ker.get.norm;
            kerNorm<1024><<<1, 1024>>>(args.size, args.arg, args.result);
            
        } else if(ker.id == EvolveL::ID) {
            
            auto args = ker.get.evolveL;
            kerEvolveL<<<args.I, args.J>>>(args.time, args.eig, args.shift, args.arg);
            
        } else if(ker.id == CopyEvolveL::ID) {
            
            auto args = ker.get.copyEvolveL;
            kerCopyEvolveL<<<args.I, args.J>>>(args.time, args.eig, args.shift, args.dest, args.source);
            
        } else if(ker.id == Trace::ID) {
            
            auto args = ker.get.trace;
            kerTrace<WarpSize><<<1, WarpSize>>>(args.dim, args.arg, args.result);
            
        } else if(ker.id == MatrixExpH::ID) {
            
            auto args = ker.get.matrixExpH;
            cudaMemsetAsync(args.arg, 0, args.dim*args.dim*sizeof(double));
            kerExpH<<<(args.dim + WarpSize - 1)/WarpSize, WarpSize>>>(args.dim, args.time, args.eig, args.shift, args.arg, args.dim + 1);
            
        } else if(ker.id == AccE::ID) {
            
            auto args = ker.get.accE;
            kerAccE<WarpSize><<<1, WarpSize>>>(args.dim, args.eig, args.arg, args.result);
            
        } else if(ker.id == Axpy::ID) {
            
            auto args = ker.get.axpy;
            cublasDaxpy(device_handle, args.size, &args.fact, args.source, 1, args.dest, 1);
        }

    }            
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

double* TrDevice::Comm::getCallBack(std::function<void(double)> callBack) {
    int index = callBack_[iStream_].size();
    callBack_[iStream_].push_back(callBack);
    return deviceCallBackBuffer_[iStream_] + index;
};

template<class K>
K& TrDevice::Comm::getKernel() {
    Kernel* kernel = hostKernelBuffer_[iStream_] + NKernel_[iStream_]++;
    kernel->id = K::ID; return *reinterpret_cast<K*>(&kernel->get);
};

void TrDevice::Comm::launch() {
    if(NKernel_[iStream_]) {
        cudaMemcpyAsync(deviceKernelBuffer_[iStream_], hostKernelBuffer_[iStream_], NKernel_[iStream_]*sizeof(Kernel), cudaMemcpyHostToDevice, stream_[iStream_]);
        kerLauncher<<<1, 1, device_, stream_[iStream_]>>>(deviceKernelBuffer_[iStream_], NKernel_[iStream_]);
        cudaMemcpyAsync(hostCallBackBuffer_[iStream_], deviceCallBackBuffer_[iStream_], callBack_[iStream_].size()*sizeof(double), cudaMemcpyDeviceToHost, stream_[iStream_]);
    } 
};

int TrDevice::Comm::is_ready(int iStream) {
    if(cudaSuccess != cudaStreamQuery(stream_[iStream])) return 0;
       
    iStream_ = iStream;
    for(std::size_t index = 0; index < callBack_[iStream_].size(); ++index)
        callBack_[iStream_][index](hostCallBackBuffer_[iStream_][index]);
    callBack_[iStream_].clear(); NKernel_[iStream_] = 0;
        
    return 1;
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TrDevice::Matrix::Matrix(double time, EigenValues const& eig) : I_(eig.dim()), J_(eig.dim()), data_(Comm::alloc(I_*J_)), exponent_(time*eig.min()) {
    auto& args = Comm::getKernel<MatrixExpH>();
    
    args.time  = time;
    args.shift = exponent_;
    args.eig   = get(eig.data());
    args.arg   = get(data_);
    args.dim   = eig.dim();
}

void TrDevice::copyEvolveL(Vector const& prop, Matrix& dest, Matrix const& source) {
    dest.I_ = source.I_; dest.J_ = source.J_; dest.exponent_ = source.exponent_ + prop.exponent(); // eigentli source.exponent_ = 0, isch aber sicherer so
    
    auto& args = Comm::getKernel<CopyEvolveL>();
    
    args.time   = prop.time();
    args.shift  = prop.exponent();
    args.eig    = get(prop.eig().data());
    args.source = get(source.data_);
    args.dest   = get(dest.data_);
    args.I      = source.I_;
    args.J      = source.J_;
}

void TrDevice::mult(Matrix& dest, Matrix const& L, Matrix const& R) {
    dest.I_ = L.I_; dest.J_ = R.J_; dest.exponent_ = L.exponent_ + R.exponent_;
    
    auto& args = Comm::getKernel<Mult>();
    
    args.A = get(L.data_);
    args.B = get(R.data_);
    args.C = get(dest.data_);
    args.M = L.I_;
    args.N = R.J_;
    args.K = L.J_;
}

void TrDevice::evolveL(Vector const& prop, Matrix& arg) {
    arg.exponent_ += prop.exponent();
    
    auto& args = Comm::getKernel<EvolveL>();
    
    args.time   = prop.time();
    args.shift  = prop.exponent();
    args.eig    = get(prop.eig().data());
    args.arg    = get(arg.data_);
    args.I      = arg.I_;
    args.J      = arg.J_;
}

void TrDevice::trace(Zahl::Zahl* Z, Zahl::Zahl* accZ, Matrix const& matrix) {
    auto& args = Comm::getKernel<Trace>(); double exponent = matrix.exponent_;
    
    args.arg    = get(matrix.data_);
    args.result = Comm::getCallBack([=](double buffer) { Zahl::Zahl temp(buffer, exponent); *Z = temp; *accZ += temp;});
    args.dim    = matrix.I_;
}

void TrDevice::accE(double* result, Zahl::Zahl fact, Matrix const& matrix, EigenValues const& eig) {
    auto& args = Comm::getKernel<AccE>(); double exponent = matrix.exponent_;
    
    args.eig    = get(eig.data());
    args.arg    = get(matrix.data_);
    args.result = Comm::getCallBack([=](double buffer) { *result += (fact*Zahl::Zahl(buffer, exponent)).toDouble();});
    args.dim    = eig.dim();
}

void TrDevice::norm(double* norm, Matrix const& matrix) {
    auto& args = Comm::getKernel<Norm>(); double exponent = matrix.exponent_;
    
    args.arg    = get(matrix.data_);
    args.result = Comm::getCallBack([=](double buffer){ *norm = std::log(buffer)/2. + exponent;});
    args.size   = matrix.I_*matrix.J_;
}


void TrDevice::axpy(Matrix& dest, Zahl::Zahl const& fact, Matrix const& source) {
    auto& args = Comm::getKernel<Axpy>();
    
    args.source    = get(source.data_);
    args.dest      = get(dest.data_);
    args.fact      = (fact*Zahl::pow(source.exponent_)).toDouble();
    args.size      = source.I_*source.J_;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__global__ void kerCublasHandle(int const flag) {
    flag ? cublasCreate(&device_handle) : cublasDestroy(device_handle);
};


void TrDevice::Comm::init(int device, int NStreams, double MiB, int bufferSize) {
    NStreams_ = NStreams; device_ = device; size_ = (MiB*(1 << 20)); bufferSize_ = bufferSize;
    
    cudaSetDevice(device_);
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
   
    ///////////////////////////////////////////////////////////////
    cudaError_t check = cudaMalloc((void**)&memory_, size_);
    if(cudaSuccess != check) {
         std::cout << cudaGetErrorString(check) << std::endl;
    }
    ///////////////////////////////////////////////////////////////
 
    stream_.resize(NStreams);
    for(int iStream = 0; iStream < NStreams; ++iStream)
        cudaStreamCreate(&stream_[iStream]);
    
    kerCublasHandle<<<1, 1>>>(1);
    cudaDeviceSynchronize();

    NKernel_.resize(NStreams, 0);
    hostKernelBuffer_.resize(NStreams);
    deviceKernelBuffer_.resize(NStreams);
    
    callBack_.resize(NStreams);
    deviceCallBackBuffer_.resize(NStreams);
    hostCallBackBuffer_.resize(NStreams);
    
    for(int iStream = 0; iStream < NStreams; ++iStream) {
        cudaMallocHost((void**)&hostKernelBuffer_[iStream], bufferSize_*sizeof(Kernel));
        deviceKernelBuffer_[iStream]   = get_aligned_memory<Kernel>(bufferSize_);
        
        cudaMallocHost((void**)&hostCallBackBuffer_[iStream], bufferSize_*sizeof(double));
        deviceCallBackBuffer_[iStream] = get_aligned_memory<double>(bufferSize_);
    }
    
    double* ptr = get_aligned_memory<double>(0);
    allocator_ = new Allocator(ptr, (size_ - pos_)/sizeof(double));
}

void TrDevice::Comm::release() {
    delete allocator_;
    
    for(int iStream = NStreams_; iStream--;) {
        cudaFreeHost(hostKernelBuffer_[iStream]);
        cudaFreeHost(hostCallBackBuffer_[iStream]);
    }
    
    cudaFree(memory_);
    
    for(int iStream = NStreams_; iStream--;)
        cudaStreamDestroy(stream_[iStream]);
    
    kerCublasHandle<<<1, 1>>>(0);
    cudaDeviceSynchronize();
}


int TrDevice::Comm::device_; int TrDevice::Comm::NStreams_;
std::size_t TrDevice::Comm::size_; std::size_t TrDevice::Comm::bufferSize_;

char* TrDevice::Comm::memory_ = 0; std::size_t TrDevice::Comm::pos_ = 0;

int TrDevice::Comm::iStream_ = 0; std::vector<cudaStream_t> TrDevice::Comm::stream_;

std::vector<int> TrDevice::Comm::NKernel_;
std::vector<Kernel*> TrDevice::Comm::hostKernelBuffer_;
std::vector<Kernel*> TrDevice::Comm::deviceKernelBuffer_;

std::vector<std::vector<std::function<void(double)> > > TrDevice::Comm::callBack_;
std::vector<double*> TrDevice::Comm::hostCallBackBuffer_;
std::vector<double*> TrDevice::Comm::deviceCallBackBuffer_;

Allocator* TrDevice::Comm::allocator_ = 0;

//-------------------------------------------------------------------------------------------------------------------------------------------------




