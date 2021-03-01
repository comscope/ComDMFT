#ifndef CTQMC_INCLUDE_UPDATES_EXPANSION_INSERTERASEPAIRCSC_H
#define CTQMC_INCLUDE_UPDATES_EXPANSION_INSERTERASEPAIRCSC_H

#include <stdexcept>
#include <iostream>
#include <vector>

#include "InsertErasePair.h"

namespace upd {
    
    namespace expansion {
        
        
        template<typename Worm>
        struct InsertCsc : Insert {
            template<typename Value>
            InsertCsc(jsx::value const& jParams, data::Data<Value> const& data) : Insert(jParams, data) {};
            
            template<typename Value>
            bool propose(double const urn, data::Data<Value> const& data, state::State<Value>& state, ut::UniformRng& urng) {
                def_ = defs_[std::upper_bound(prob_.begin(), prob_.end(), urn*prob_.back()) - prob_.begin()];
                
                auto const& opsL = state.expansion()[def_.flavorL];
                auto const& opsR = state.expansion()[def_.flavorR];
                
                ut::KeyType keyHigh;
                if(opsL.size() && opsR.size()) {
                    std::size_t const index = urng()*(opsL.size() + opsR.size());
                    ut::KeyType const keyLow = index < opsL.size() ? opsL[index] : opsR[index - opsL.size()];
                    
                    auto it = std::upper_bound(opsL.begin(), opsL.end(), keyLow);
                    keyHigh = it != opsL.end() ? *it : *opsL.begin() + ut::KeyMax;
                    
                    it = std::upper_bound(opsR.begin(), opsR.end(), keyLow);
                    keyHigh = std::min(it != opsR.end() ? *it : *opsR.begin() + ut::KeyMax, keyHigh);
                    
                    keyDiff_ = keyHigh - keyLow; if(keyHigh > ut::KeyMax) keyHigh -= ut::KeyMax;
                } else
                    keyHigh = keyDiff_ = ut::KeyMax;
                
                def_.keyL = ut::cyclic(keyHigh - urng()*keyDiff_);
                def_.keyR = ut::cyclic(keyHigh - urng()*keyDiff_);

                return true;
            };
            
            Worm const& target(Worm const& worm) const {
                return worm;
            };
            
            template<typename Value>
            double ratio(data::Data<Value> const& data, state::State<Value>& state) const {
                int const N = state.expansion()[def_.flavorL].size() + state.expansion()[def_.flavorR].size();
                double const timeDiff = (ut::beta()*keyDiff_)/ut::KeyMax;
                return timeDiff*timeDiff*(N ? N/(N + 2.) : 1.);
            };
            
        private:
            ut::KeyType keyDiff_;
        };
        
        
        template<typename Worm>
        struct EraseCsc : Erase {
            template<typename Value>
            EraseCsc(jsx::value const& jParams, data::Data<Value> const& data) : Erase(jParams, data) {};
            
            template<typename Value>
            bool propose(double const urn, data::Data<Value> const& data, state::State<Value>& state, ut::UniformRng& urng) {
                def_ = defs_[std::upper_bound(prob_.begin(), prob_.end(), urn*prob_.back()) - prob_.begin()];
                
                auto const& opsL = state.expansion()[def_.flavorL];
                auto const& opsR = state.expansion()[def_.flavorR];
                
                if(!(opsL.size() && opsR.size())) return false;
                
                if(opsL.size() + opsR.size() > 2) {
                    std::size_t index = urng()*(opsL.size() + opsR.size());
                    ut::KeyType const keyLow = index < opsL.size() ? opsL[index] : opsR[index - opsL.size()];
                    
                    auto itL = std::upper_bound(opsL.begin(), opsL.end(), keyLow); ut::KeyType shiftL = 0;
                    if(itL != opsL.end()) def_.keyL = *itL; else { def_.keyL = *(itL = opsL.begin()); shiftL = ut::KeyMax;};
                    
                    auto itR = std::upper_bound(opsR.begin(), opsR.end(), keyLow); ut::KeyType shiftR = 0;
                    if(itR != opsR.end()) def_.keyR = *itR; else { def_.keyR = *(itR = opsR.begin()); shiftR = ut::KeyMax;};
                    
                    ut::KeyType keyHigh = ++itL != opsL.end() ? *itL + shiftL : *opsL.begin() + ut::KeyMax;
                    keyHigh = std::min(++itR != opsR.end() ? *itR + shiftR : *opsR.begin() + ut::KeyMax, keyHigh);
                    
                    if(std::max(def_.keyL + shiftL, def_.keyR + shiftR) < keyHigh) keyDiff_ = keyHigh - keyLow; else return 0;
                } else {
                    def_.keyL = *opsL.begin();
                    def_.keyR = *opsR.begin();
                    keyDiff_ = ut::KeyMax;
                }

                return true;
            };
            
            Worm const& target(Worm const& worm) const {
                return worm;
            };
            
            template<typename Value>
            double ratio(data::Data<Value> const& data, state::State<Value>& state) const {
                int const N = state.expansion()[def_.flavorL].size() + state.expansion()[def_.flavorR].size();
                double const timeDiff = (ut::beta()*keyDiff_)/ut::KeyMax;
                return (N > 2 ? N/(N - 2.) : 1.)/(timeDiff*timeDiff);
            };
            
        private:
            ut::KeyType keyDiff_;
        };

        
    }
    
}

#endif