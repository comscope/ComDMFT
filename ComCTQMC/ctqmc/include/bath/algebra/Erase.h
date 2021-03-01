#ifndef CTQMC_INCLUDE_BATH_ALGEBRA_ERASEPAIR_H
#define CTQMC_INCLUDE_BATH_ALGEBRA_ERASEPAIR_H

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
    
    struct Erase {
        ut::KeyType keyL;
        ut::KeyType keyR;
    };
    
    template<typename Value>
    struct Update<Value, Erase> : itf::Update<Value> {
        
        TypeId type() {
            return get_type_id<Erase>();
        };
        
        void add(Erase upd, Bath<Value> const& bath, Hyb<Value> const& hyb) {
            if(guard_) throw std::runtime_error("bath::Update<Erase>: use multi-erase, not multiple erases");
            
            itL_ = bath.posL_.find(upd.keyL); itR_ = bath.posR_.find(upd.keyR); guard_ = true;
        };
        
        double ratio(Bath<Value> const& bath, Hyb<Value> const& hyb) {
            return std::abs(val_ = bath.B_.at(itR_->second, itL_->second));
        };
        
        int accept(Bath<Value>& bath, Hyb<Value> const& hyb) {
            int const N = bath.opsL_.size(); int const newN = N - 1;
            int const posL = itL_->second; int const posR = itR_->second;
            
            int sign = 1;  Matrix<Value> temp(newN);
            
            if(newN) {
                int const inc = 1;
        
                if(posL != newN) {
                    swap(&N, bath.B_.data(0, newN), &inc, bath.B_.data(0, posL), &inc); sign *= -1;
                    bath.posL_[bath.opsL_.back().key()] = posL; bath.opsL_[posL] = bath.opsL_.back();
                }
                if(posR != newN) {
                    swap(&N, bath.B_.data(newN, 0), &N, bath.B_.data(posR, 0), &N); sign *= -1;
                    bath.posR_[bath.opsR_.back().key()] = posR; bath.opsR_[posR] = bath.opsR_.back();
                }
                
                for(int n = 0; n < newN; ++n)
                    copy(&newN, bath.B_.data(0, n), &inc, temp.data(0, n), &inc);
                
                Value const fact = -1./val_;
                geru(&newN, &newN, &fact, bath.B_.data(0, newN), &inc, bath.B_.data(newN, 0), &N, temp.data(), &newN);
            }
            
            bath.B_ = std::move(temp);
            
            bath.posL_.erase(itL_); bath.opsL_.pop_back();
            bath.posR_.erase(itR_); bath.opsR_.pop_back();
            
            bath.det_ *= val_*static_cast<double>(sign); return sign;
        };
        
        void reject(Bath<Value>& bath, Hyb<Value> const& hyb) {
        };
        
    private:
        bool guard_ = false; Value val_;
        
        std::map<ut::KeyType, int>::const_iterator itL_;
        std::map<ut::KeyType, int>::const_iterator itR_;
    };
}

#endif
