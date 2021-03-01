#include <cstring>
#include <type_traits>
#include <sstream>

#include "Algebra.h"
#include "Variant.h"

#include "../include/Errchk.h"
#include "../include/Allocator.h"

#include "../../../include/BlasLapack.h"
#include "../../../include/mpi/Basic.h"

#ifdef HAVE_CUBLAS

#include <cublasLt.h>

__device__ double deviceDZero = .0;
__device__ double deviceDOne  = 1.;
__device__ cuComplex deviceCZero = {.0,.0};
__device__ cuComplex deviceCOne  = {1.,.0};
__device__ cublasLtHandle_t deviceHandle;

void cublasHandle(int const flag) {
    flag ? cublasLtCreate(&deviceHandle) : cublasLtDestroy(deviceHandle);
};

#else

#include <cutlass/gemm/gemm.h>
#include <cutlass/gemm/dgemm_traits.h>
/*
 //TODO: Try half precision -- this won't compile on summit
 #include <cutlass/gemm/volta884_complex_gemm_traits.h>
using cutlass_planar_complex_traits = typename cutlass::gemm::Volta884ComplexGemmTraits<
        cutlass::MatrixLayout::kColumnMajor, cutlass::MatrixTransform::kNone,
        cutlass::MatrixLayout::kColumnMajor, cutlass::MatrixTransform::kNone,
        cutlass::Shape<8, 64, 128>, cutlass::Shape<4, 32, 32>, //
        double,double,double,
        2>;
*/
#endif


using namespace imp;
using namespace device;


constexpr int WarpSize = 32;  Allocator* alloc = nullptr;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int imp::pci_id_size() {
    return 16;
}

// @#$!#$@!# cudaGetDevicePciBusId is not working properly on SummitDev ......
void get_pci_id(char* pci_id, int deviceId) {
    std::stringstream stream;
    cudaDeviceProp deviceProperties; cudaErrchk(cudaGetDeviceProperties(&deviceProperties, deviceId));
    stream << std::hex << deviceProperties.pciDomainID << ":" << deviceProperties.pciBusID << ":" << deviceProperties.pciDeviceID;
    std::string str = stream.str(); std::copy(str.begin(), str.end(), pci_id);
}


int imp::get_pci_ids(std::vector<char>& pciIds) {
    int deviceCount=0;
    cudaErrchk(cudaGetDeviceCount(&deviceCount));
        
    pciIds.resize(deviceCount*pci_id_size(), '\0');
    for(int id = 0; id < deviceCount; ++id) get_pci_id(&pciIds[id*pci_id_size()], id);
    
    return deviceCount;
}

void imp::init_device(int const deviceId, std::size_t processesPerDevice) {
    //int deviceId = mpi::rank_on_node(); //cudaErrchk(cudaDeviceGetByPCIBusId(&deviceId, pciId.data()));
    cudaErrchk(cudaSetDevice(deviceId));
    
    cudaDeviceProp deviceProperties; cudaErrchk(cudaGetDeviceProperties(&deviceProperties, deviceId));
    if(deviceProperties.computeMode != cudaComputeModeExclusive && deviceProperties.computeMode != cudaComputeModeExclusiveProcess)
        throw std::runtime_error("Please set GPU compute mode to \"cudaComputeModeExclusive\" or \"cudaComputeModeExclusiveProcess\"");
    if(deviceProperties.warpSize != WarpSize)
        throw std::runtime_error("Please set WarpSize in AlgebraDevice.cu to " + std::to_string(deviceProperties.warpSize));
    
    cudaErrchk(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));
    
#ifdef HAVE_CUBLAS
    
    CublasHandle(1);
    
#endif
    
    alloc = new Allocator((0.8/processesPerDevice)*deviceProperties.totalGlobalMem);
    
    cudaErrchk(cudaDeviceSynchronize());
}

void imp::release_device() {
    if(!alloc->sanity_check()) throw std::runtime_error("Memory leak !");
    
    delete alloc;
    alloc = nullptr;
    
#ifdef HAVE_CUBLAS
    
    CublasHandle(0);
    
#endif
    
    cudaErrchk(cudaDeviceSynchronize());
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
template <typename Value>
imp::Matrix<Device, Value>::Matrix(int size) : data_(alloc->get<cuda_value_t<Value>>(size * cuda_value<Value>::size)) {
}

template <typename Value>
imp::Matrix<Device, Value>::Matrix(Matrix<Device, Value>::Identity const& identity) : I_(identity.dim), J_(identity.dim), data_(alloc->get<cuda_value_t<Value>>(cuda_value<Value>::size*I_*J_)), exponent_(.0) {
    
    cuda_value_t<Value>* temp = new cuda_value_t<Value>[cuda_value<Value>::size*I_*J_];
    
    std::memset(temp, 0, cuda_value<Value>::size*I_*J_*sizeof(cuda_value_t<Value>));
    
    for(int i = 0; i < identity.dim; ++i) temp[i*(identity.dim + 1)] = 1.; //kack memset isch das allgemein für double's ?
    
    cudaErrchk(cudaMemcpy(data_.ptr(), temp, cuda_value<Value>::size*I_*J_*sizeof(cuda_value_t<Value>), cudaMemcpyHostToDevice));
    
    delete[] temp;
}

template <typename Value>
imp::Matrix<Device, Value>::Matrix(Matrix<Device, Value>::Zero const& zero) : I_(zero.dim), J_(zero.dim), data_(alloc->get<cuda_value_t<Value>>(cuda_value<Value>::size*I_*J_)), exponent_(.0) {
    
    cuda_value_t<Value>* temp = new cuda_value_t<Value>[cuda_value<Value>::size*I_*J_];
    
    std::memset(temp, 0, cuda_value<Value>::size*I_*J_*sizeof(cuda_value_t<Value>));
    
    cudaErrchk(cudaMemcpy(data_.ptr(), temp, cuda_value<Value>::size*I_*J_*sizeof(cuda_value_t<Value>), cudaMemcpyHostToDevice));
    
    delete[] temp;
}

namespace imp{

    template <>
    Matrix<Device, double>::Matrix(int I, int J, io::Matrix<double> const& mat) : I_(I), J_(J), data_(alloc->get<cuda_value_t<double>>(I_*J_)), exponent_(.0) {
        double* temp = new double[I_*J_];
        
        for(int i = 0; i < I; ++i)
            for(int j = 0; j < J; ++j)
                temp[j + J*i] = mat(i, j);
        
        cudaErrchk(cudaMemcpy(data_.ptr(), temp, I_*J_*sizeof(double), cudaMemcpyHostToDevice));
        
        delete[] temp;
    }
        
    template <>
    Matrix<Device, ut::complex>::Matrix(int I, int J, io::Matrix<ut::complex> const& mat) : I_(I), J_(J), data_(alloc->get<cuda_value_t<ut::complex>>(cuda_value<ut::complex>::size*I_*J_)), exponent_(.0) {
        cuda_value_t<ut::complex>* temp = new cuda_value_t<ut::complex>[cuda_value<ut::complex>::size*I_*J_];
        
        //planar complex structure
        for(int i = 0; i < I; ++i)
            for(int j = 0; j < J; ++j){
                temp[j + J*i] = mat(i, j).real();
                temp[j + J*i + I*J] = mat(i, j).imag();
            }
        
        cudaErrchk(cudaMemcpy(data_.ptr(), temp, cuda_value<ut::complex>::size*I_*J_*sizeof(cuda_value_t<ut::complex>), cudaMemcpyHostToDevice));
        
        delete[] temp;
    }

}
    
template <typename Value>
imp::Matrix<Device, Value>::~Matrix() {
    alloc->free(data_);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template <>
void imp::add<double>(double* dest, double fact, Matrix<Device, double> const& source) {
    int const N = source.I()*source.J(); int const one = 1;
    cuda_value_t<double>* temp = new cuda_value_t<double>[N];
    
    cudaErrchk(cudaMemcpy(temp, source.data().ptr(), N*sizeof(cuda_value_t<double>), cudaMemcpyDeviceToHost));
    
    daxpy_(&N, &fact, temp, &one, dest, &one);
    
    delete[] temp;
}

//TODO: do this with cudaMemcpy2d so that the cuda planar complex -> ut::complex
template <>
void imp::add<ut::complex>(ut::complex* dest, ut::complex fact, Matrix<Device, ut::complex> const& source) {
    int const N = source.I()*source.J(); int const one = 1;
    
    cuda_value_t<ut::complex>* planar_complex_temp = new cuda_value_t<ut::complex>[cuda_value<ut::complex>::size*N];
    ut::complex* ut_complex_temp = new ut::complex[N];
    
    cudaErrchk(cudaMemcpy(planar_complex_temp, source.data().ptr(), cuda_value<ut::complex>::size*N*sizeof(cuda_value_t<ut::complex>), cudaMemcpyDeviceToHost));
    
    for (int i=0; i<N; i++)
        ut_complex_temp[i] = ut::complex(planar_complex_temp[i], planar_complex_temp[i + N]);
    
    zaxpy_(&N, &fact, ut_complex_temp, &one, dest, &one);

    delete[] planar_complex_temp;
    delete[] ut_complex_temp;
    
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//-------------------------------------------------------------------------------------------------------------------------------------------------


template <typename Value>
struct CopyEvolveL {
    double time;
    double shift;
    double const* energies;
    cuda_value_t<Value> const* source;
    cuda_value_t<Value>* dest;
    int I;
    int J;
    int size;
};


__global__ void kerCopyEvolveL(CopyEvolveL<double> args) {
    int const i = blockIdx.x; int const j = threadIdx.x;

    args.dest[j + blockDim.x*i] = exp(args.time*args.energies[i] - args.shift)*args.source[j + blockDim.x*i];
};

__global__ void kerCopyEvolveL(CopyEvolveL<ut::complex> args) {
    int const i = blockIdx.x; int const j = threadIdx.x;

    args.dest[j + blockDim.x*i] = exp(args.time*args.energies[i] - args.shift)*args.source[j + blockDim.x*i];
    args.dest[j + blockDim.x*i + args.size] = exp(args.time*args.energies[i] - args.shift)*args.source[j + blockDim.x*i + args.size];
    
};


template <typename Value>
void imp::copyEvolveL(Matrix<Device, Value>& dest, Vector<Device> const& prop, Matrix<Device, Value> const& source, itf::Batcher<Value>& batcher) {
    dest.I() = source.I(); dest.J() = source.J(); dest.exponent() = source.exponent() + prop.exponent(); // eigentli source.exponent() = 0, isch aber sicherer so
    
    auto& args = imp::get<Device>(batcher).template get_kernel<CopyEvolveL<Value>>();
    
    args.time     = prop.time();
    args.shift    = prop.exponent();
    args.energies = prop.energies().data().ptr();
    args.source   = source.data().ptr();
    args.dest     = dest.data().ptr();
    args.I        = source.I();
    args.J        = source.J();
    args.size     = source.I()*source.J(); //only needed for complex
    
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template <typename Value>
struct Mult {
    cuda_value_t<Value> const* A_real;
    cuda_value_t<Value> const* B_real;
    cuda_value_t<Value>* C_real;
    int M;
    int N;
    int K;
    
    //pointers to second half of A/B/C matrices, i.e., imaginary part
    cuda_value_t<Value> const* A_imag;
    cuda_value_t<Value> const* B_imag;
    cuda_value_t<Value>* C_imag;
    
#ifdef HAVE_CUBLAS
    
    cublasStatus_t status;
    cublasLtMatmulDesc_t matmulDesc;
    cublasLtMatrixLayout_t Adesc, Bdesc, Cdesc;
    
    int64_t AplaneOffset;
    int64_t BplaneOffset;
    int64_t CplaneOffset;
    
#endif
    
};

#ifndef HAVE_CUBLAS

template <typename KernelClass>
__global__ void cutlass_kernel(typename KernelClass::Params const& params)
{
    extern __shared__ int GemmSharedStorageBase[];
    
    typename KernelClass::SharedStorage *shared_storage =
    reinterpret_cast<typename KernelClass::SharedStorage *>(GemmSharedStorageBase);
    
    KernelClass gemm(params, *shared_storage);
    
    gemm.multiply_add();
}

template<typename BlockShape, typename ThreadShape>
__device__ void cutlass_gemm(Mult<double> const& args, Byte*& memory)
{
    typedef cutlass::gemm::DgemmTraits<
    cutlass::MatrixLayout::kColumnMajor,   // layout of A matrix
    cutlass::MatrixLayout::kColumnMajor,
    BlockShape,
    cutlass::gemm::LinearScaling<double>,
    ThreadShape
    >
    Traits;
    
    typedef typename Traits::Params Params;
    typedef typename Traits::KernelClass KernelClass;

    memory = reinterpret_cast<Byte*>(reinterpret_cast<unsigned long long>(memory + (alignof(Params) - 1)) & -alignof(Params));

    // Params needs to be trivially destructible ... do not see how to test this in compile time (there is no equivalent to std::is_trivially_destructible in thrust as far as I can see
    Params& params = *new(memory) Params();  memory += sizeof(Params);

    params.initialize(
                      args.M,
                      args.N,
                      args.K,
                      1.,
                      args.A_real,
                      args.M,
                      args.B_real,
                      args.K,
                      .0,
                      args.C_real,
                      args.M,
                      args.C_real,
                      args.M
                      );
    
    cutlass_kernel<KernelClass><<< params.grid, params.block, sizeof(typename KernelClass::SharedStorage)>>>(params);
};


template<typename BlockShape, typename ThreadShape>
__device__ void cutlass_gemm(Mult<ut::complex> const& args, Byte*& memory)
{
    typedef cutlass::gemm::DgemmTraits<
       cutlass::MatrixLayout::kColumnMajor,   // layout of A matrix
       cutlass::MatrixLayout::kColumnMajor,
       BlockShape,
       cutlass::gemm::LinearScaling<double>,
       ThreadShape
       > Traits;
    typedef typename Traits::Params Params;
    typedef typename Traits::KernelClass KernelClass;

    memory = reinterpret_cast<Byte*>(reinterpret_cast<unsigned long long>(memory + (alignof(Params) - 1)) & -alignof(Params));

    // Params needs to be trivially destructible ... do not see how to test this in compile time (there is no equivalent to std::is_trivially_destructible in thrust as far as I can see
    Params& params_rr = *new(memory) Params();  memory += sizeof(Params);
    Params& params_ii = *new(memory) Params();  memory += sizeof(Params);
    Params& params_ri = *new(memory) Params();  memory += sizeof(Params);
    Params& params_ir = *new(memory) Params();  memory += sizeof(Params);


    params_rr.initialize(
                      args.M,
                      args.N,
                      args.K,
                      1.,
                      args.A_real,
                      args.M,
                      args.B_real,
                      args.K,
                      .0,
                      args.C_real,
                      args.M,
                      args.C_real,
                      args.M
                      );
    
    params_ii.initialize(
                      args.M,
                      args.N,
                      args.K,
                      -1.,
                      args.A_imag,
                      args.M,
                      args.B_imag,
                      args.K,
                      1.,
                      args.C_real,
                      args.M,
                      args.C_real,
                      args.M
                      );
    
    params_ri.initialize(
                      args.M,
                      args.N,
                      args.K,
                      1.,
                      args.A_imag,
                      args.M,
                      args.B_real,
                      args.K,
                      .0,
                      args.C_imag,
                      args.M,
                      args.C_imag,
                      args.M
                      );
    
    params_ir.initialize(
                      args.M,
                      args.N,
                      args.K,
                      1.,
                      args.A_real,
                      args.M,
                      args.B_imag,
                      args.K,
                      1.,
                      args.C_imag,
                      args.M,
                      args.C_imag,
                      args.M
                      );
    
    
    cutlass_kernel<KernelClass><<< params_rr.grid, params_rr.block, sizeof(typename KernelClass::SharedStorage)>>>(params_rr);
    cutlass_kernel<KernelClass><<< params_ri.grid, params_ri.block, sizeof(typename KernelClass::SharedStorage)>>>(params_ri);
    
    __syncthreads(); //needed?
    
    cutlass_kernel<KernelClass><<< params_ii.grid, params_ii.block, sizeof(typename KernelClass::SharedStorage)>>>(params_ii);
    cutlass_kernel<KernelClass><<< params_ir.grid, params_ir.block, sizeof(typename KernelClass::SharedStorage)>>>(params_ir);
    
};

#endif

template <typename Value>
void imp::mult(Matrix<Device, Value>& dest, Matrix<Device, Value> const& L, Matrix<Device, Value> const& R, itf::Batcher<Value>& batcher) {
    dest.I() = L.I(); dest.J() = R.J(); dest.exponent() = L.exponent() + R.exponent();
    
    auto& args = imp::get<Device>(batcher).template get_kernel<Mult<Value>>();
    
    args.A_real = R.data().ptr();
    args.B_real = L.data().ptr();
    args.C_real = dest.data().ptr();
    args.M = R.J();
    args.N = L.I();
    args.K = L.J();
    
    args.A_imag = R.data().ptr()+R.I()*R.J();
    args.B_imag = L.data().ptr()+L.I()*L.J();
    args.C_imag = dest.data().ptr()+dest.I()*dest.J();
    
#ifdef HAVE_CUBLAS
    
    /*
     These functions can't be called from kernel...
    args.AplaneOffset = (args.A_imag - args.A_real) * sizeof(args.A_real[0]);
    args.BplaneOffset = (args.B_imag - args.B_real) * sizeof(args.B_real[0]);
    args.CplaneOffset = (args.C_imag - args.C_real) * sizeof(args.C_real[0]);
    
    args.status = cublasLtMatmulDescCreate(&args.matmulDesc, CUDA_C_32F);
    
    args.status = cublasLtMatrixLayoutCreate(&args.Adesc, CUDA_C_16F, args.M, args.K, args.M); //m,k,lda
    if(args.status != CUBLAS_STATUS_SUCCESS) throw std::runtime_error("CUBLAS GEMM: Error creating A matrix layout");
    args.status = cublasLtMatrixLayoutSetAttribute(args.Adesc, CUBLASLT_MATRIX_LAYOUT_PLANE_OFFSET, &args.AplaneOffset, sizeof(args.AplaneOffset));
    if(args.status != CUBLAS_STATUS_SUCCESS) throw std::runtime_error("CUBLAS GEMM: Error creating A matrix atributes");

    args.status = cublasLtMatrixLayoutCreate(&args.Bdesc, CUDA_C_16F, args.K, args.N, args.K); //k,n,lda
    if(args.status != CUBLAS_STATUS_SUCCESS) throw std::runtime_error("CUBLAS GEMM: Error creating B matrix layout");
    args.status = cublasLtMatrixLayoutSetAttribute(args.Bdesc, CUBLASLT_MATRIX_LAYOUT_PLANE_OFFSET, &args.BplaneOffset, sizeof(args.BplaneOffset));
    if(args.status != CUBLAS_STATUS_SUCCESS) throw std::runtime_error("CUBLAS GEMM: Error creating B matrix atributes");

    args.status = cublasLtMatrixLayoutCreate(&args.Cdesc, CUDA_C_16F, args.M, args.N, args.M); //m,n,ldc
    if(args.status != CUBLAS_STATUS_SUCCESS) throw std::runtime_error("CUBLAS GEMM: Error creating C matrix layout");
    args.status = cublasLtMatrixLayoutSetAttribute(args.Cdesc, CUBLASLT_MATRIX_LAYOUT_PLANE_OFFSET, &args.CplaneOffset, sizeof(args.CplaneOffset));
    if(args.status != CUBLAS_STATUS_SUCCESS) throw std::runtime_error("CUBLAS GEMM: Error creating D matrix layout");
    assert(args.status == CUBLAS_STATUS_SUCCESS);
     */
    
    
#endif
    
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template <typename Value>
struct EvolveL {
    double time;
    double shift;
    double const* energies;
    cuda_value_t<Value>* arg;
    int I;
    int J;
    int size;
};

__global__ void kerEvolveL(EvolveL<double> args) {
    int const i = blockIdx.x; int const j = threadIdx.x;
    
    args.arg[j + blockDim.x*i] *= exp(args.time*args.energies[i] - args.shift);
};

__global__ void kerEvolveL(EvolveL<ut::complex> args) {
    int const i = blockIdx.x; int const j = threadIdx.x;
    
    args.arg[j + blockDim.x*i] *= exp(args.time*args.energies[i] - args.shift);
    args.arg[j + blockDim.x*i + args.size] *= exp(args.time*args.energies[i] - args.shift);
    
};


template <typename Value>
void imp::evolveL(Vector<Device> const& prop, Matrix<Device, Value>& arg, itf::Batcher<Value>& batcher) {
    arg.exponent() += prop.exponent();
    
    auto& args = imp::get<Device>(batcher).template get_kernel<EvolveL<Value>>();
    
    args.time     = prop.time();
    args.shift    = prop.exponent();
    args.energies = prop.energies().data().ptr();
    args.arg      = arg.data().ptr();
    args.I        = arg.I();
    args.J        = arg.J();
    args.size     = arg.I()*arg.J();
    
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

#if __CUDACC_VER_MAJOR__ >= 9

template <typename T>
__device__ __forceinline__ void reduceWarp(int const tid, T* data, T* result) {
    T temp;
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

template <typename T>
__device__ __forceinline__ void reduceWarp(int const tid, T volatile* data, T* result) {
    data[tid] += data[tid + 16];
    data[tid] += data[tid + 8];
    data[tid] += data[tid + 4];
    data[tid] += data[tid + 2];
    data[tid] += data[tid + 1];
    if(tid == 0) *result = *data;
};

#endif

template <int Size>
__device__ __forceinline__ void reduce(int const tid, double* data, double* result) {
    if(tid < Size/2) data[tid] += data[tid + Size/2];
    __syncthreads();
    reduce<Size/2>(tid, data, result);
};

template <>
__device__ __forceinline__ void reduce<WarpSize>(int const tid, double* data, double* result) {
    if(tid < WarpSize) reduceWarp<double>(tid, data, result);
};


template <int Size>
__device__ __forceinline__ void reduce(int const tid, float* data, float* result) {
    if(tid < Size/2) data[tid] += data[tid + Size/2];
    __syncthreads();
    reduce<Size/2>(tid, data, result);
};

template <>
__device__ __forceinline__ void reduce<WarpSize>(int const tid, float* data, float* result) {
    if(tid < WarpSize) reduceWarp<float>(tid, data, result);
};

template <int Size>
__device__ __forceinline__ void reduce(int const tid, thrust::complex<double>* data, thrust::complex<double>* result) {
    if(tid < Size/2) data[tid] += data[tid + Size/2];
    __syncthreads();
    reduce<Size/2>(tid, data, result);
};

template <>
__device__ __forceinline__ void reduce<WarpSize>(int const tid, thrust::complex<double>* data, thrust::complex<double>* result) {
    if(tid < WarpSize) reduceWarp<thrust::complex<double>>(tid, data, result);
};

template <int Size>
__device__ __forceinline__ void reduce(int const tid, thrust::complex<float>* data, thrust::complex<float>* result) {
    if(tid < Size/2) data[tid] += data[tid + Size/2];
    __syncthreads();
    reduce<Size/2>(tid, data, result);
};

template <>
__device__ __forceinline__ void reduce<WarpSize>(int const tid, thrust::complex<float>* data, thrust::complex<float>* result) {
    if(tid < WarpSize) reduceWarp<thrust::complex<float>>(tid, data, result);
};



//-------------------------------------------------------------------------------------------------------------------------------------------------

template <typename Value>
struct Trace {
    cuda_value_t<Value> const* arg;
    cuda_value_scalar<Value>* result;
    int dim;
    int size; // only needed for complex
};

template <int BlockDim>
__global__ void kerTrace(Trace<double> args) {
    __shared__ cuda_value_scalar<double> cache[BlockDim + 16];   // I do not want some threads in the reduceWarp to read stuff outside the cache ... nobody of the nVidia freaks seems to care about this (and probabely they are right) but I do not see why
    cache[threadIdx.x] = .0;
    int i = threadIdx.x;
    
    while(i < args.dim) {
        cache[threadIdx.x] += args.arg[(args.dim + 1)*i];
        
        i += BlockDim;
    }
    
    __syncthreads();
    
    reduce<BlockDim>(threadIdx.x, cache, args.result);
};

template <int BlockDim>
__global__ void kerTrace(Trace<ut::complex> args) {
    __shared__ cuda_value_scalar<ut::complex> cache[BlockDim + 16];
    cache[threadIdx.x] = .0;
    int i = threadIdx.x;
    
    while(i < args.dim) {
        
        cache[threadIdx.x] += cuda_value_scalar<ut::complex>(
                                                            args.arg[(args.dim + 1)*i],
                                                            args.arg[(args.dim + 1)*i + args.size]
                                                            );

        i += BlockDim;
    }
    
    __syncthreads();
    
    reduce<BlockDim>(threadIdx.x, cache, args.result);
};


template <typename Value>
void imp::trace(ut::Zahl<Value>* Z, ut::Zahl<Value>* accZ, Matrix<Device, Value> const& matrix, itf::Batcher<Value>& batcher) {
    auto& args = imp::get<Device>(batcher).template get_kernel<Trace<Value>>(); double exponent = matrix.exponent();
    
    args.arg    = matrix.data().ptr();
    args.result = imp::get<Device>(batcher).get_callback([=](cuda_value_scalar<Value> buffer) { ut::Zahl<Value> temp(buffer, exponent); if(Z) *Z = temp; if(accZ) *accZ += temp;});
    args.dim    = matrix.I();
    args.size   = matrix.I()*matrix.J();
    
}


//-------------------------------------------------------------------------------------------------------------------------------------------------

template <typename Value>
struct TraceAtB {
    cuda_value_t<Value> const* At;
    cuda_value_t<Value> const* B;
    cuda_value_scalar<Value>* result;
    int size;
};

template <int BlockDim>
__global__ void kerTraceAtB(TraceAtB<double> args) {
    __shared__ cuda_value_scalar<double> cache[BlockDim + 16];
    cache[threadIdx.x] = .0;
    int i = threadIdx.x;
    
    while(i < args.size) {
        cache[threadIdx.x] += args.At[i]*args.B[i];
        
        i += BlockDim;
    }
    
    __syncthreads();
    
    reduce<BlockDim>(threadIdx.x, cache, args.result);
};

template <int BlockDim>
__global__ void kerTraceAtB(TraceAtB<ut::complex> args) {
    __shared__ cuda_value_scalar<ut::complex> cache[BlockDim + 16];
    cache[threadIdx.x] = .0;
    int i = threadIdx.x;
    
    while(i < args.size) {
        

        cache[threadIdx.x] += cuda_value_scalar<ut::complex>(
                                                             args.At[i            ] * args.B[i            ] //real x real
                                                            -args.At[i + args.size] * args.B[i + args.size] //imag x imag
                                                            ,args.At[i + args.size] * args.B[i            ] //imag x real
                                                            +args.At[i            ] * args.B[i + args.size] //real x imag
                                                                 );
            

        
        i += BlockDim;
    }
    
    __syncthreads();
    
    reduce<BlockDim>(threadIdx.x, cache, args.result);
};



template <typename Value>
void imp::traceAtB(ut::Zahl<Value>* Z, ut::Zahl<Value>* accZ, Matrix<Device, Value> const& At, Matrix<Device, Value> const& B, itf::Batcher<Value>& batcher) {
    auto& args = imp::get<Device>(batcher).template get_kernel<TraceAtB<Value>>(); double exponent = At.exponent() + B.exponent();
    
    args.At     = At.data().ptr();
    args.B      = B.data().ptr();
    args.result = imp::get<Device>(batcher).get_callback([=](cuda_value_scalar<Value> buffer) { ut::Zahl<Value> temp(buffer, exponent); if(Z) *Z = temp; if(accZ) *accZ += temp;});
    args.size   = At.I()*At.J();

}


//-------------------------------------------------------------------------------------------------------------------------------------------------

template <typename Value>
struct Norm {
    cuda_value_t<Value> const* arg;
    cuda_value_scalar<Value>* result;
    int size;
};

template<int BlockDim>
__global__ void kerNorm(Norm<double> args) {
    __shared__ cuda_value_scalar<double> cache[BlockDim + 16];
    cache[threadIdx.x] = .0;
    int i = threadIdx.x;
    
    while(i < args.size) {
        cache[threadIdx.x] += args.arg[i]*args.arg[i]; 
        
        i += BlockDim;
    }
    
    __syncthreads();
    
    reduce<BlockDim>(threadIdx.x, cache, args.result);
};

template<int BlockDim>
__global__ void kerNorm(Norm<ut::complex> args) {
    __shared__ cuda_value_scalar<ut::complex> cache[BlockDim + 16];
    cache[threadIdx.x] = .0;
    int i = threadIdx.x;
    
    while(i < args.size) {

        cache[threadIdx.x] += args.arg[i]*args.arg[i] + args.arg[i + args.size]*args.arg[i + args.size];

        i += BlockDim;
    }
    
    __syncthreads();
    
    reduce<BlockDim>(threadIdx.x, cache, args.result);
};

namespace imp{
    template <>
    void norm(double* norm, Matrix<Device, double> const& matrix, itf::Batcher<double>& batcher) {
        auto& args = imp::get<Device>(batcher).template get_kernel<Norm<double>>(); double exponent = matrix.exponent();
        
        args.arg    = matrix.data().ptr();
        args.result = imp::get<Device>(batcher).get_callback([=](cuda_value_scalar<double> buffer) { *norm = std::log(buffer)/2. + exponent;});
        args.size   = matrix.I()*matrix.J();
        
    }

    template <>
    void norm(double* norm, Matrix<Device, ut::complex> const& matrix, itf::Batcher<ut::complex>& batcher) {
        auto& args = imp::get<Device>(batcher).template get_kernel<Norm<ut::complex>>(); double exponent = matrix.exponent();
        
        args.arg    = matrix.data().ptr();
        args.result = imp::get<Device>(batcher).get_callback([=](cuda_value_scalar<ut::complex> buffer) { *norm = std::log(thrust::abs(buffer))/2. + exponent;});
        args.size   = matrix.I()*matrix.J();
        
    }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template <typename Value>
void imp::density_matrix(Matrix<Device, Value>& dest, Matrix<Device, Value> const& B, Vector<Device> const& prop, Matrix<Device, Value> const& A, Energies<Device> const& energies, itf::Batcher<Value>& batcher)
{
    throw std::runtime_error("imp::density_matrix: not implemented !");
};

//-------------------------------------------------------------------------------------------------------------------------------------------------

//CUDA doesn't like non-pointer thrust::complex because it has a copy constructor.
//So we have to use the C version
template <typename Value>
struct Add {
    cuda_value_t<Value> const* source;
    cuda_value_t<Value>* dest;
    cuda_value_cscalar<Value> fact;
    int size;
};


__global__ void kerAdd(Add<double> args)
{
    int const index = blockDim.x*blockIdx.x + threadIdx.x;
    
    if(index < args.size) args.dest[index] += args.fact*args.source[index];
};

__global__ void kerAdd(Add<ut::complex> args)
{
    int const index = blockDim.x*blockIdx.x + threadIdx.x;
    
    //could give this to 2x or 4x threads. Probably a small part of GPU time though
    if(index < args.size){
    
        args.dest[index]            += cuCreal(args.fact)*args.source[index];           //real x real
        args.dest[index]            -= cuCimag(args.fact)*args.source[index+args.size]; //imag * imag
        args.dest[index+args.size]  += cuCimag(args.fact)*args.source[index];           //imag * real
        args.dest[index+args.size]  += cuCreal(args.fact)*args.source[index+args.size]; //real * imag
        
    }
    
};

//need namespace because of specialization
namespace imp{
template <>
void add(Matrix<Device, double>& dest, ut::Zahl<double> const& fact, Matrix<Device, double> const& source, itf::Batcher<double>& batcher)
{
    auto& args = imp::get<Device>(batcher).template get_kernel<Add<double>>();
    
    args.source    = source.data().ptr();
    args.dest      = dest.data().ptr();
    args.fact      = (fact*ut::exp(source.exponent())).get();
    args.size      = source.I()*source.J();
    
}

template <>
void add(Matrix<Device, ut::complex>& dest, ut::Zahl< ut::complex> const& fact, Matrix<Device,  ut::complex> const& source, itf::Batcher< ut::complex>& batcher)
{
    auto& args = imp::get<Device>(batcher).template get_kernel<Add<ut::complex>>();
    
    auto const temp = (fact*ut::Zahl<ut::complex>(1., source.exponent())).get();
    
    args.source    = source.data().ptr();
    args.dest      = dest.data().ptr();
    args.fact      = make_cuDoubleComplex(temp.real(), temp.imag());
    args.size      = source.I()*source.J();
    
}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename Value>
using KerArgs = variant<CopyEvolveL<Value>, Mult<Value>, EvolveL<Value>, Trace<Value>, TraceAtB<Value>, Norm<Value>, Add<Value>>;

//alignof(cuda_value_t<Value>) better than 16?
template <typename Value>
struct alignas(16) imp::Kernel {
    KerArgs<Value> args;
    int id;
};



template <typename Value>
__global__ void kerLauncher(Kernel<Value>* kernel, int const N, Byte* memory)
{
    

    for(int n = 0; n < N; ++n) {
        Kernel<Value>&  ker = kernel[n];
        
        if(ker.id == device::index<Mult<Value>, KerArgs<Value>>::value) {

            auto& args = get_device<Mult<Value>>(ker.args);
            
#ifdef HAVE_CUBLAS
            
            assert(("CUBLAS Not implemented", false)) // not implemented

#else
            
            cutlass_gemm<cutlass::Shape<8, 64, 128>, cutlass::Shape<8, 8, 8>>(args, memory);
            
#endif
            
        } else if(ker.id == device::index<Norm<Value>, KerArgs<Value>>::value) {

            auto& args = get_device<Norm<Value>>(ker.args);
            kerNorm<1024><<<1, 1024>>>(args);

        } else if(ker.id == device::index<EvolveL<Value>, KerArgs<Value>>::value) {

            auto& args = get_device<EvolveL<Value>>(ker.args);
            kerEvolveL<<<args.I, args.J>>>(args);

        } else if(ker.id == device::index<CopyEvolveL<Value>, KerArgs<Value>>::value) {

            auto& args = get_device<CopyEvolveL<Value>>(ker.args);
            kerCopyEvolveL<<<args.I, args.J>>>(args);

        } else if(ker.id == device::index<Trace<Value>, KerArgs<Value>>::value) {

            auto& args = get_device<Trace<Value>>(ker.args);
            kerTrace<WarpSize><<<1, WarpSize>>>(args);

        } else if(ker.id == device::index<TraceAtB<Value>, KerArgs<Value>>::value) {

            auto& args = get_device<TraceAtB<Value>>(ker.args);
            kerTraceAtB<1024><<<1, 1024>>>(args);

        } else if(ker.id == device::index<Add<Value>, KerArgs<Value>>::value) {
            
            auto& args = get_device<Add<Value>>(ker.args);
            kerAdd<<<(args.size + 256 - 1)/256, 256>>>(args);
            
        }
        
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename Value>
imp::Batcher<Device, Value>::Batcher(std::size_t size) :
size_(size),
numberOfKernels_(0),
phase_(Phase::record),
deviceKernelBuffer_(alloc->get<Kernel<Value>>(size_)),
deviceCallBackBuffer_(alloc->get<cuda_value_scalar<Value>>(size_)),
memory_(alloc->get<Byte>(
#ifndef HAVE_CUBLAS
                         sizeof(typename cutlass::gemm::DgemmTraits<
                                cutlass::MatrixLayout::kColumnMajor,
                                cutlass::MatrixLayout::kColumnMajor>::Params) // size of the dgemm instructions
                         *size_*cuda_value<Value>::size // size of the matrices being multiplied -- why -- x2 right??
                         *cuda_value<Value>::size*cuda_value<Value>::size // have to store 4 params for planar complex dgemm
#else
                         8
#endif
                         ))
{
    cudaErrchk(cudaStreamCreate(&stream_));
    cudaErrchk(cudaMallocHost(reinterpret_cast<void**>(&hostKernelBuffer_), cuda_value<Value>::size*size_*sizeof(Kernel<Value>)));
    cudaErrchk(cudaMallocHost(reinterpret_cast<void**>(&hostCallBackBuffer_), cuda_value<Value>::size*size_*sizeof(cuda_value_scalar<Value>)));
}

template <typename Value>
cuda_value_scalar<Value>* imp::Batcher<Device, Value>::get_callback(std::function<void(cuda_value_scalar<Value>)> callBack) {
    if(phase_ != Phase::record) throw std::runtime_error("imp::Batcher::get_callback");
    
    int index = callBack_.size();  callBack_.push_back(callBack);
    return deviceCallBackBuffer_.ptr() + index;
};

template <typename Value>
template <typename K>
K& imp::Batcher<Device, Value>::get_kernel() {
    if(phase_ != Phase::record) throw std::runtime_error("imp::Batcher::get_kernel");
    
    Kernel<Value>& ker = hostKernelBuffer_[numberOfKernels_++];
    ker.id = device::index<K, KerArgs<Value>>::value;
    return get_host<K>(ker.args);
};

template <typename Value>
void imp::Batcher<Device, Value>::launch() {
    if(phase_ != Phase::record) throw std::runtime_error("imp::Batcher::launch");
    
    if(numberOfKernels_) {
        cudaErrchk(cudaMemcpyAsync(deviceKernelBuffer_.ptr(), hostKernelBuffer_, numberOfKernels_*sizeof(Kernel<Value>), cudaMemcpyHostToDevice, stream_));
        kerLauncher<Value><<<1, 1, 0, stream_>>>(deviceKernelBuffer_.ptr(), numberOfKernels_, memory_.ptr());
        
        numberOfKernels_ = 0; phase_ = Phase::execute;
    }
};

template <typename Value>
int imp::Batcher<Device, Value>::is_ready() {
    if(phase_ == Phase::execute) {
        cudaError_t quest = cudaStreamQuery(stream_);
        
        if(quest == cudaErrorNotReady) return 0;
        
        cudaErrchk(quest);

        if(callBack_.size()) {
            cudaErrchk(cudaMemcpyAsync(hostCallBackBuffer_, deviceCallBackBuffer_.ptr(), callBack_.size()*sizeof(cuda_value_scalar<Value>), cudaMemcpyDeviceToHost, stream_));
            phase_ = Phase::finalize; return 0;
        }
    }
    
    if(phase_ == Phase::finalize) {
        cudaError_t quest = cudaStreamQuery(stream_);
        
        if(quest == cudaErrorNotReady) return 0;
        
        cudaErrchk(quest);
        
        for(std::size_t index = 0; index < callBack_.size(); ++index) callBack_[index](hostCallBackBuffer_[index]);
        callBack_.clear();
    }
    
    phase_ = Phase::record; return 1;
};

template <typename Value>
imp::Batcher<Device, Value>::~Batcher() {
    alloc->free(memory_);
    
    alloc->free(deviceCallBackBuffer_);
    cudaErrchk(cudaFreeHost(hostCallBackBuffer_));
    
    alloc->free(deviceKernelBuffer_);
    cudaErrchk(cudaFreeHost(hostKernelBuffer_));
    
    cudaErrchk(cudaStreamDestroy(stream_));
};


//explicit instantiations (double)
namespace imp{
    template struct Matrix<Device,double>;
    template struct Batcher<Device,double>;
    template void copyEvolveL(Matrix<Device, double>& dest, Vector<Device> const& prop, Matrix<Device, double> const& source, itf::Batcher<double>& batcher);
    template void mult(Matrix<Device, double>& dest, Matrix<Device, double> const& L, Matrix<Device, double> const& R, itf::Batcher<double>& batcher);
    template void evolveL(Vector<Device> const& prop, Matrix<Device, double>& arg, itf::Batcher<double>& batcher);
    template void trace(ut::Zahl<double>* Z, ut::Zahl<double>* accZ, Matrix<Device, double> const& matrix, itf::Batcher<double>& batcher);
    template void traceAtB(ut::Zahl<double>* Z, ut::Zahl<double>* accZ, Matrix<Device, double> const& At, Matrix<Device, double> const& B, itf::Batcher<double>& batcher);
    template void norm(double* norm, Matrix<Device, double> const& matrix, itf::Batcher<double>& batcher);
    template void add(Matrix<Device, double>& dest, ut::Zahl<double> const& fact, Matrix<Device, double> const& source, itf::Batcher<double>& batcher);
    //template void add(double* dest, double fact, Matrix<Device, double> const& source);
    template void density_matrix(Matrix<Device, double>& dest, Matrix<Device, double> const& B, Vector<Device> const& prop, Matrix<Device, double> const& A, Energies<Device> const& energies, itf::Batcher<double>& batcher);


    //explicit instantiations (complex)
    template struct Matrix<Device,ut::complex>;
    template struct Batcher<Device,ut::complex>;
    template void copyEvolveL(Matrix<Device, ut::complex>& dest, Vector<Device> const& prop, Matrix<Device, ut::complex> const& source, itf::Batcher<ut::complex>& batcher);
    template void mult(Matrix<Device, ut::complex>& dest, Matrix<Device, ut::complex> const& L, Matrix<Device, ut::complex> const& R, itf::Batcher<ut::complex>& batcher);
    template void evolveL(Vector<Device> const& prop, Matrix<Device, ut::complex>& arg, itf::Batcher<ut::complex>& batcher);
    template void trace(ut::Zahl<ut::complex>* Z, ut::Zahl<ut::complex>* accZ, Matrix<Device, ut::complex> const& matrix, itf::Batcher<ut::complex>& batcher);
    template void traceAtB(ut::Zahl<ut::complex>* Z, ut::Zahl<ut::complex>* accZ, Matrix<Device, ut::complex> const& At, Matrix<Device, ut::complex> const& B, itf::Batcher<ut::complex>& batcher);
    template void norm(double* norm, Matrix<Device, ut::complex> const& matrix, itf::Batcher<ut::complex>& batcher);
    template void add(Matrix<Device, ut::complex>& dest, ut::Zahl<ut::complex> const& fact, Matrix<Device, ut::complex> const& source, itf::Batcher<ut::complex>& batcher);
    //template void add(ut::complex* dest, ut::complex fact, Matrix<Device, ut::complex> const& source);
    template void density_matrix(Matrix<Device, ut::complex>& dest, Matrix<Device, ut::complex> const& B, Vector<Device> const& prop, Matrix<Device, ut::complex> const& A, Energies<Device> const& energies, itf::Batcher<ut::complex>& batcher);
}
