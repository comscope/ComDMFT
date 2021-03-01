#ifndef CTQMC_INCLUDE_BATH_ALGEBRA_EXCHANGE_BATH_H
#define CTQMC_INCLUDE_BATH_ALGEBRA_EXCHANGE_BATH_H

#include <cmath>
#include <stdexcept>
#include <vector>
#include <map>
#include <algorithm>

#include "../Bath.h"
#include "../Hyb.h"
#include "../../../../include/BlasLapack.h"


namespace bath {
    
    struct Exchange {
        ut::KeyType keyOld, keyNew;
        int flavorNew;
    };
    
    template<typename Value>
    struct Update<Value, Exchange> : itf::Update<Value> {

        TypeId type() {
            return get_type_id<Exchange>();
        };
       
        void add(Exchange upd, Bath<Value> const& bath, Hyb<Value> const& hyb) {
            int const N = bath.opsL_.size(); std::vector<Value> vec(N);

            if(hyb.isR(upd.flavorNew)) {                                                                    // u = Delta, v = e
                int const pos = bath.posR_.find(upd.keyOld)->second, flavorOld = bath.opsR_[pos].flavor();  auto const& ops = bath.opsL_;

                for(int n = 0; n < N; ++n)
                    vec[n] = hyb(ops[n].flavor(), upd.flavorNew, ops[n].key() - upd.keyNew) - hyb(ops[n].flavor(), flavorOld, ops[n].key() - upd.keyOld);
                
                for(auto& list : list_)
                    if(!list.isR) list.vec[pos] = vec[list.pos] = .5*(hyb(list.flavorNew, upd.flavorNew, list.keyNew - upd.keyNew) - hyb(ops[list.pos].flavor(), flavorOld, ops[list.pos].key() - upd.keyOld));
                
                list_.push_back({true, pos, upd.keyOld, upd.keyNew, upd.flavorNew, std::move(vec)});
            } else {                                                                                    // u = e, v = Delta
                int const pos = bath.posL_.find(upd.keyOld)->second, flavorOld = bath.opsL_[pos].flavor();  auto const& ops = bath.opsR_;
                
                for(int n = 0; n < N; ++n)
                    vec[n] = hyb(upd.flavorNew, ops[n].flavor(), upd.keyNew - ops[n].key()) - hyb(flavorOld, ops[n].flavor(), upd.keyOld - ops[n].key());
                
                for(auto& list : list_)
                    if(list.isR) list.vec[pos] = vec[list.pos] = .5*(hyb(upd.flavorNew, list.flavorNew, upd.keyNew - list.keyNew) - hyb(flavorOld, ops[list.pos].flavor(), upd.keyOld - ops[list.pos].key()));
                
                list_.push_back({false, pos, upd.keyOld, upd.keyNew, upd.flavorNew, std::move(vec)});
            }
        };
        
        double ratio(Bath<Value> const& bath, Hyb<Value> const& hyb) {
            int const inc = 1, N = bath.opsL_.size(), size = list_.size(); auto const& B = bath.B_;
            
            if(size == 1) {
                auto const& list = list_[0];
                return std::abs(detRatio_ = list.isR ? 1. + dotu(&N, B.data(list.pos, 0), &N, list.vec.data(), &inc) : 1. +  dotu(&N, list.vec.data(), &inc, B.data(0, list.pos), &inc));
            }
            
            std::vector<Value> toInvert(size*size); auto it = toInvert.data();
            for(auto const& u_j : list_)
                if(u_j.isR) {
                    for(auto const& v_i : list_)
                        if(v_i.isR) {                                                                  // v_i = e,      u_j = Delta
                            *it++ = dotu(&N, B.data(v_i.pos, 0), &N, u_j.vec.data(), &inc);
                        } else {                                                                       // v_i = Delta,  u_j = Delta
                            char const no = 'n'; Value const zero = .0, one = 1.; std::vector<Value> Bu(N);
                            gemv(&no, &N, &N, &one, B.data(), &N, u_j.vec.data(), &inc, &zero, Bu.data(), &inc);
                            *it++ = dotu(&N, v_i.vec.data(), &inc, Bu.data(), &inc);
                        }
                } else
                    for(auto const& v_i : list_)                                                       // v_i = e,      u_j = e     :     v_i = Delta,  u_j = e
                        *it++ = v_i.isR ? B.at(v_i.pos, u_j.pos) : dotu(&N, v_i.vec.data(), &inc, B.data(0, u_j.pos), &inc);
            
            for(int i = 0; i < size; ++i) toInvert[i*size + i] += 1.;
            
            inv_.resize(size*size, .0); for(int i = 0; i < size; ++i) inv_[i*size + i] = 1.;
            
            int ipiv[size]; int info;
            gesv(&size, &size, toInvert.data(), &size, ipiv, inv_.data(), &size, &info);
            
            detRatio_ = 1.;
            for(int i = 0; i < size; ++i)
                detRatio_ *= (ipiv[i] != i + 1 ? -toInvert[i*size + i] : toInvert[i*size + i]);
            
            return std::abs(detRatio_);
        };
        
        int accept(Bath<Value>& bath, Hyb<Value> const& hyb) {
            int const inc = 1, N = bath.opsL_.size(), size = list_.size(); auto& B = bath.B_;
            char const no = 'n', yes = 't'; Value const zero = .0, one = 1.;
            
            std::vector<Value> BU(N*size), BtV(N*size);
            for(int k = 0; k < size; ++k) {
                auto const& list = list_[k];
                if(list.isR) {                                                                   // u = Delta, v = e
                    gemv(&no, &N, &N, &one, B.data(), &N, list.vec.data(), &inc, &zero, BU.data() + N*k, &inc);
                    copy(&N, B.data(list.pos, 0), &N, BtV.data() + N*k, &inc);
                    
                    bath.posR_.erase(list.keyOld); bath.posR_[list.keyNew] = list.pos;
                    bath.opsR_[list.pos] = Operator<Value>(list.keyNew, list.flavorNew);
                } else {                                                                         // u = e,     v = Delta
                    copy(&N, B.data(0, list.pos), &inc, BU.data() + N*k, &inc);
                    gemv(&yes, &N, &N, &one, B.data(), &N, list.vec.data(), &inc, &zero, BtV.data() + N*k, &inc);
                    
                    bath.posL_.erase(list.keyOld); bath.posL_[list.keyNew] = list.pos;
                    bath.opsL_[list.pos] = Operator<Value>(list.keyNew, list.flavorNew);
                }
            }
            
            if(size == 1) {
                Value const fact = -1./detRatio_;
                geru(&N, &N, &fact, BU.data(), &inc, BtV.data(), &inc, B.data(), &N);
            } else {
                Value const fact = -1.;  std::vector<Value> temp(N*size);
                gemm(&no, &no, &N, &size, &size, &one, BU.data(), &N, inv_.data(), &size, &zero, temp.data(), &N);
                gemm(&no, &yes, &N, &N, &size, &fact, temp.data(), &N, BtV.data(), &N, &one, B.data(), &N);
            }
            
            bath.det_ *= detRatio_;  return 1;
        };
        
        void reject(Bath<Value>& bath, Hyb<Value> const& hyb) {
        };
        
    private:
        struct Entry {
            bool isR;
            int pos;
            ut::KeyType keyOld, keyNew, flavorNew;
            std::vector<Value> vec;
        };
                                                           
        Value detRatio_;
        std::vector<Value> inv_;
        std::vector<Entry> list_;
        
    };
    
}

#endif
