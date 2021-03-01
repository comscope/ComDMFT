#ifndef CTQMC_DEVICE_INCLUDE_ALLOCATOR_H
#define CTQMC_DEVICE_INCLUDE_ALLOCATOR_H

#include <stdexcept>
#include <list>
#include <set>
#include <cassert>
#include <limits>
#include <iostream>
#include <vector>
#include <algorithm>

#include "Errchk.h"

#include "../../include/Utilities.h"


namespace device {
    
    using Byte = unsigned char;
    
    class Allocator {
        enum struct State { free, occupied };
        
        struct DataImpl {
            Byte* ptr;
            std::size_t bytes;
            State state;
            
            DataImpl() = delete;
            DataImpl(Byte* ptr, std::size_t bytes, State state) : ptr(ptr), bytes(bytes), state(state) {};
            ~DataImpl() = default;
        };
        
        typedef std::list<DataImpl>::iterator PtrImpl;
        
        struct LessPtrImpl {
            bool operator()(PtrImpl const& lhs, PtrImpl const& rhs) {
                if(lhs->bytes < rhs->bytes) return true;
                if(lhs->bytes > rhs->bytes) return false;
                return lhs->ptr < rhs->ptr;
            }
        };
        
    public:
        template<typename T>
        struct Data {
            Data() = delete;    // delete some other constructors as well, assignments ???
            Data(PtrImpl ptrImpl) :
            ptrImpl_(ptrImpl), ptr_(reinterpret_cast<T*>(reinterpret_cast<std::uintptr_t>(ptrImpl_->ptr + (alignof(T) - 1)) & -alignof(T))) {
            };
            
            T* ptr() { return ptr_;};
            T const* ptr() const { return ptr_;};
        private:
            PtrImpl ptrImpl_;
            T* const ptr_;
            
            friend class Allocator;
        };
        
        Allocator() = delete;
        Allocator(std::size_t bytes) : bytes_(bytes), memory_(nullptr)  {
            cudaErrchk(cudaMalloc(reinterpret_cast<void**>(&memory_), bytes_));

            //assert that memory aligned to alignof(double), just for the sake of cleaness
            
            all_.emplace_back(nullptr, 0, State::occupied);
            free_.insert(all_.emplace(all_.end(), memory_, bytes_, State::free));
            all_.emplace_back(nullptr, 0, State::occupied);
            
            dummy_ = all_.emplace(all_.end(), nullptr, 0, State::free);  // ??????
        };
        Allocator(Allocator const&) = delete;
        Allocator(Allocator&&) = delete;
        Allocator& operator=(Allocator const&) = delete;
        Allocator& operator=(Allocator&&) = delete;
        
        template<typename T>
        Data<T> get(std::size_t size) {
            if(!size) throw std::runtime_error("Allocation of 0 bytes !"); // do we need this ???

            dummy_->bytes = sizeof(T)*size + alignof(T) - 1;

            std::set<PtrImpl>::iterator it = free_.lower_bound(dummy_);
            if(it == free_.end()) throw ut::out_of_memory();
            
            PtrImpl temp = *it; free_.erase(it);
            
            if(temp->bytes == dummy_->bytes) {
                temp->state = State::occupied;
                return temp;
            }
            
            PtrImpl cut = all_.emplace(temp, temp->ptr, dummy_->bytes, State::occupied);
            
            temp->bytes -= dummy_->bytes;
            temp->ptr += dummy_->bytes;
            
            free_.insert(temp);

            return cut;
        }

        template<typename T>
        void free(Data<T> data) {
            if(data.ptrImpl_ == all_.begin()) return;

            PtrImpl next = data.ptrImpl_; ++next;
            if(next->state == State::free) {
                free_.erase(next);
                data.ptrImpl_->bytes += next->bytes;
                all_.erase(next);
            }
            
            PtrImpl prev = data.ptrImpl_; --prev;
            if(prev->state == State::free) {
                free_.erase(prev);
                data.ptrImpl_->ptr = prev->ptr;
                data.ptrImpl_->bytes += prev->bytes;
                all_.erase(prev);
            }
            
            data.ptrImpl_->state = State::free;
            free_.insert(data.ptrImpl_);
        }

        bool sanity_check() {
            return free_.size() == 1 && (*free_.begin())->ptr == memory_ && (*free_.begin())->bytes == bytes_; 
        }
        
        ~Allocator() {
            cudaErrchk(cudaFree(memory_));
        };
    private:
        std::size_t const bytes_;
        Byte* memory_;

        PtrImpl dummy_;
        
        std::list<DataImpl> all_;
        std::set<PtrImpl, LessPtrImpl> free_;
    };
}


#endif
