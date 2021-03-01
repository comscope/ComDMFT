#ifndef CTQMC_INCLUDE_BATH_ALGEBRA_INSERTPAIR_H
#define CTQMC_INCLUDE_BATH_ALGEBRA_INSERTPAIR_H

#include <cmath>
#include <stdexcept>
#include <vector>
#include <map>
#include <algorithm>

#include "../Bath.h"
#include "../Hyb.h"
#include "../../../../include/BlasLapack.h"

//Vielleicht insert und erase noch verändern wie in CTBI hürütütütüt pluschter pluschter

namespace bath {
    
    struct Insert{
        ut::KeyType keyL; int flavorL;
        ut::KeyType keyR; int flavorR;
    };
    
    template<typename Value>
    struct Update<Value, Insert> : itf::Update<Value> {
        
        TypeId type() {
            return get_type_id<Insert>();
        };
        
        void add(Insert upd,  Bath<Value> const& bath, Hyb<Value> const& hyb) {
            if(guard_) throw std::runtime_error("bath::Update<Insert>: use multi-erase, not multiple erases");
            
            int const N = bath.opsL_.size();
            
            val_ = hyb(upd.flavorL, upd.flavorR, upd.keyL - upd.keyR);
            
            if(N) {
                Bv_.resize(N); vec_.resize(N);
                for(int n = 0; n < N; ++n) vec_[n] = hyb(bath.opsL_[n].flavor(), upd.flavorR, bath.opsL_[n].key() - upd.keyR);
                
                char const no = 'n';
                int const inc = 1;
                Value const zero = .0;
                Value const one = 1.;
                gemv(&no, &N, &N, &one, bath.B_.data(), &N, vec_.data(), &inc, &zero, Bv_.data(), &inc);
                
                for(int n = 0; n < N; ++n) vec_[n] = hyb(upd.flavorL, bath.opsR_[n].flavor(), upd.keyL - bath.opsR_[n].key());
                val_ -= dotu(&N, vec_.data(), &inc, Bv_.data(), &inc);
            }
            
            upd_ = upd;  guard_ = true;
        };
        
        double ratio(Bath<Value> const& bath, Hyb<Value> const& hyb) {
            return std::abs(val_);
        };

        int accept(Bath<Value>& bath, Hyb<Value> const& hyb) {
            int const N = bath.opsL_.size();
            int const newN = N + 1;
            
            bath.posL_[upd_.keyL] = bath.opsL_.size(); bath.opsL_.push_back(Operator<Value>(upd_.keyL, upd_.flavorL));
            bath.posR_[upd_.keyR] = bath.opsR_.size(); bath.opsR_.push_back(Operator<Value>(upd_.keyR, upd_.flavorR));
            
            Matrix<Value> temp(newN);
            
            Value fact = 1./val_;
            temp.at(N, N) = fact;
            
            if(N) {
                std::vector<Value> hBTilde(N);
                
                char const yes = 't';
                int const inc = 1;
                Value const zero = .0;
                Value const one = 1.;
                gemv(&yes, &N, &N, &fact, bath.B_.data(), &N, vec_.data(), &inc, &zero, hBTilde.data(), &inc);
                geru(&N, &N, &one, Bv_.data(), &inc, hBTilde.data(), &inc, bath.B_.data(), &N);
                
                for(int n = 0; n < N; ++n)
                    copy(&N, bath.B_.data(0, n), &inc, temp.data(0, n), &inc);
                
                fact = -1./val_;
                scal(&N, &fact, Bv_.data(), &inc);
                copy(&N, Bv_.data(), &inc, temp.data(0, N), &inc);
                
                Value const minus = -1.;
                scal(&N, &minus, hBTilde.data(), &inc);
                copy(&N, hBTilde.data(), &inc, temp.data(N, 0), &newN);
            }
            
            bath.B_ = std::move(temp);  bath.det_ *= val_; return 1;
        };
        
        void reject(Bath<Value>& bath, Hyb<Value> const& hyb) {
        };

    private:
        bool guard_ = false; Insert upd_;
        
        Value val_;
        std::vector<Value> Bv_;
        std::vector<Value> vec_;
    };
    
}

#endif
