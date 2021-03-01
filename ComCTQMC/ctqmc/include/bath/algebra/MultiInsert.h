#ifndef CTQMC_INCLUDE_BATH_ALGEBRA_MULIT_INSERT_H
#define CTQMC_INCLUDE_BATH_ALGEBRA_MULIT_INSERT_H

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
    
    struct MultiInsert{
        std::vector<Insert> VecInsert;
    };
    
    template<typename Value>
    struct Update<Value, MultiInsert> : itf::Update<Value> {
        
        TypeId type() {
            return get_type_id<MultiInsert>();
        };
        
        void add(MultiInsert updates,  Bath<Value> const& bath, Hyb<Value> const& hyb) {
            if(guard_) throw std::runtime_error("bath::Update<Insert>: use multi-insert, not multiple insert");
            
            updates_ = updates;  guard_ = true;
        };
        
        double ratio(Bath<Value> const& bath, Hyb<Value> const& hyb) {
            
            target_ = bath;
            
            for (auto const& upd : updates_.VecInsert){
                target_.insertL(upd.keyL, upd.flavorL);
                target_.insertR(upd.keyR, upd.flavorR);
            }
            
            target_.clean(hyb); //better way than recomputing the whole matrix?
            
            return std::abs(target_.det_/bath.det_);
        };

        int accept(Bath<Value>& bath, Hyb<Value> const& hyb) {
            
            bath=target_;
            
            return 1;
        };
        
        void reject(Bath<Value>& bath, Hyb<Value> const& hyb) {
        };

    private:
        bool guard_ = false;
        MultiInsert updates_;
        Bath<Value> target_;
    };
    
}

#endif
