#ifndef CTQMC_INCLUDE_BATH_ALGEBRA_MULTI_ERASE_H
#define CTQMC_INCLUDE_BATH_ALGEBRA_MULTI_ERASE_H

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
    
    struct MultiErase {
        std::vector<Erase> VecErase;
    };
    
    template<typename Value>
    struct Update<Value, MultiErase> : itf::Update<Value> {
        
        TypeId type() {
            return get_type_id<MultiErase>();
        };
        
        void add(MultiErase updates, bath::Bath<Value>& bath, Hyb<Value> const& hyb) {
            if(guard_) throw std::runtime_error("bath::Update<Erase>: use multi-erase, not multiple erases");
            
            updates_ = updates;
            guard_ = true;
        };
        
        double ratio(Bath<Value> const& bath, Hyb<Value> const& hyb) {
            
            sign = 1;
            target_ = bath;
            for (auto const& upd : updates_.VecErase){
                sign*=target_.eraseL(upd.keyL);
                sign*=target_.eraseR(upd.keyR);
            }
            
            target_.clean(hyb); //better way than recomputing the whole matrix?
            
            return std::abs(target_.det_/bath.det_);
        };
        
        int accept(Bath<Value>& bath, Hyb<Value> const& hyb) {
            
            bath = std::move(target_);
            
            return sign; //sign flip
        };
        
        void reject(Bath<Value>& bath, Hyb<Value> const& hyb) {
        };
        
    private:
        bool guard_ = false;
        MultiErase updates_;
        Bath<Value> target_;
        int sign;
    };
}

#endif
