#ifndef TRACEALGEBRADEVICE
#define TRACEALGEBRADEVICE

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

#include "../include/Utilities.h"

namespace TrDevice {
    struct Allocator {
        enum struct State { free, occupied };
        
        struct Data {
            struct LessPtr {
                bool operator()(std::list<Data>::iterator const& lhs, std::list<Data>::iterator const& rhs) {
                    if(lhs->size_ < rhs->size_) return true;
                    if(lhs->size_ > rhs->size_) return false;
                    return lhs->data_ < rhs->data_;
                }
            };
            
            Data(double* data, std::size_t size, State state) : data_(data), size_(size), state_(state) {};
            ~Data() {};
        private:
            double* data_;
            std::size_t size_;
            State state_;
            
            friend struct Allocator;
            friend double* get(std::list<Data>::iterator);
        };
        
        typedef std::list<Data>::iterator DataPtr;
        
        
        
        Allocator(double* memory, std::size_t size) : memory_(memory), size_(size) {
            all_.emplace_back(static_cast<double*>(0), 0, State::occupied);
            free_.insert(all_.emplace(all_.end(), memory_, size_, State::free));
            all_.emplace_back(static_cast<double*>(0), 0, State::occupied);
            
            dummy_ = all_.emplace(all_.end(), static_cast<double*>(0), 0, State::free);
        };

        
        DataPtr alloc(std::size_t size) {
            if(!size) throw std::runtime_error("Allocation of 0 bytes !");
            
            dummy_->size_ = size;
            
            std::set<DataPtr>::iterator it = free_.lower_bound(dummy_);
            if(it == free_.end()) throw Ut::out_of_memory();
            
            DataPtr temp = *it; free_.erase(it);
            
            if(temp->size_ == size) {
                temp->state_ = State::occupied;
                return temp;
            }
            
            DataPtr cut = all_.emplace(temp, temp->data_, size, State::occupied);
            
            temp->size_ -= size;
            temp->data_ += size;
            
            free_.insert(temp);

            return cut;
        }
        
        void free(std::list<Data>::iterator it) {
            DataPtr next = it; ++next;
            if(next->state_ == State::free) {
                free_.erase(next);
                it->size_ += next->size_;
                all_.erase(next);
            }
            
            DataPtr prev = it; --prev;
            if(prev->state_ == State::free) {
                free_.erase(prev);
                it->data_ = prev->data_;
                it->size_ += prev->size_;
                all_.erase(prev);
            }
            
            it->state_ = State::free;
            free_.insert(it);
        }
        
        
        std::size_t remaining() {
            if(free_.size() != 1) throw std::runtime_error("Error during initialisation !");
            
            return (*free_.begin())->size_;
        }
        
        ~Allocator() {
            if(free_.size() != 1)
                throw std::runtime_error("Memory Leak");
            
            if((*free_.begin())->data_ != memory_ || (*free_.begin())->size_ != size_)
                throw std::runtime_error("Memory Leak");
        }
    private:
        double* const memory_;
        std::size_t const size_;
        
        DataPtr dummy_;
        
        std::list<Data> all_;
        std::set<DataPtr, Data::LessPtr> free_;
    };
    
    typedef Allocator::DataPtr data_ptr; inline double* get(data_ptr ptr) { return ptr->data_;};
    
    
    struct Kernel;
    
    struct Comm {
        template<class T>
        static T* get_aligned_memory(std::size_t size) {
            if(256%alignof(T)) throw std::runtime_error("TrDevice::Comm::get_aligned_memory: can not ensure proper alignement !");
            pos_ = ((pos_ + alignof(T) - 1)/alignof(T))*alignof(T);
            T* ptr = reinterpret_cast<T*>(memory_ + pos_); pos_ += size*sizeof(T);
            if(pos_ > size_) throw std::runtime_error("TrDevice::Comm::get_aligned_memory: out of memory !");
            return ptr;
        }
        
        static void init(int device, int NStreams, double MiB, int bufferSize);
        static void release();
        
        static Allocator::DataPtr alloc(std::size_t size) { return allocator_->alloc(size);};
        static void free(Allocator::DataPtr ptr) { allocator_->free(ptr);};

        static double* getCallBack(std::function<void(double)> callBack);
        template<class T> static T& getKernel();
        
        static int is_ready(int);
        static void launch();
      
        static int iStream() { return iStream_;}; 
    private:
        static int device_; static int NStreams_;
        static std::size_t size_; static std::size_t bufferSize_;
        
        static char* memory_; static std::size_t pos_;
        
        static int iStream_; static std::vector<cudaStream_t> stream_;
        
        static std::vector<int> NKernel_;
        static std::vector<Kernel*> hostKernelBuffer_;
        static std::vector<Kernel*> deviceKernelBuffer_;
        
        static std::vector<std::vector<std::function<void(double)> > > callBack_;
        static std::vector<double*> hostCallBackBuffer_;
        static std::vector<double*> deviceCallBackBuffer_;
        
        static Allocator* allocator_;
    };
}


#endif
