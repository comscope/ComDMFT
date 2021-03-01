#ifndef CTQMC_INCLUDE_UPDATES_EXPANSION_INSERTERASEPAIRSTD_H
#define CTQMC_INCLUDE_UPDATES_EXPANSION_INSERTERASEPAIRSTD_H


#include "InsertErasePair.h"


namespace upd {
    
    namespace expansion {
        
        
        template<typename Worm>
        struct InsertStd : Insert {
            template<typename Value>
            InsertStd(jsx::value const& jParams, data::Data<Value> const& data) : Insert(jParams, data) {};
            
            template<typename Value>
            bool propose(double const urn, data::Data<Value> const& data, state::State<Value>& state, ut::UniformRng& urng) {
                def_ = defs_[std::upper_bound(prob_.begin(), prob_.end(), urn*prob_.back()) - prob_.begin()];
                
                def_.keyL = ut::KeyMax*urng();
                def_.keyR = ut::KeyMax*urng();
                
                return true;
            };
            
            Worm const& target(Worm const& worm) const {
                return worm;
            };
            
            template<typename Value>
            double ratio(data::Data<Value> const& data, state::State<Value>& state) const {
                return ut::beta()*ut::beta()/((state.expansion()[def_.flavorL].size() + 1)*(state.expansion()[def_.flavorR].size() + 1));
            };
        };
        
        
        template<typename Worm>
        struct EraseStd : Erase {
            template<typename Value>
            EraseStd(jsx::value const& jParams, data::Data<Value> const& data) : Erase(jParams, data) {};
            
            template<typename Value>
            bool propose(double const urn, data::Data<Value> const& data, state::State<Value>& state, ut::UniformRng& urng) {
                def_ = defs_[std::upper_bound(prob_.begin(), prob_.end(), urn*prob_.back()) - prob_.begin()];
                
                auto const& opsL = state.expansion()[def_.flavorL];
                auto const& opsR = state.expansion()[def_.flavorR];
                
                if(!(opsL.size() && opsR.size())) return false;
                
                def_.keyL = opsL[urng()*opsL.size()];
                def_.keyR = opsR[urng()*opsR.size()];
                
                return true;
            };
            
            Worm const& target(Worm const& worm) const {
                return worm;
            };
            
            template<typename Value>
            double ratio(data::Data<Value> const& data, state::State<Value>& state) const {
                return (state.expansion()[def_.flavorL].size()*state.expansion()[def_.flavorR].size())/(ut::beta()*ut::beta());
            };
            
        };
        
    }
    
}

#endif