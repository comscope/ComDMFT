#ifndef CTQMC_INCLUDE_UPDATES_EXPANSION_INSERTERASEQUADSTD_H
#define CTQMC_INCLUDE_UPDATES_EXPANSION_INSERTERASEQUADSTD_H


#include "InsertEraseQuad.h"


namespace upd {
    
    namespace expansion {
        
        
        template<typename Worm>
        struct QuadInsertStd : QuadInsert {
            template<typename Value>
            QuadInsertStd(jsx::value const& jParams, data::Data<Value> const& data) : QuadInsert(jParams, data) {};
            
            template<typename Value>
            bool propose(double const urn, data::Data<Value> const& data, state::State<Value>& state, ut::UniformRng& urng) {
                def_ = defs_[std::upper_bound(prob_.begin(), prob_.end(), urn*prob_.back()) - prob_.begin()];
                
                def_.keyL1 = ut::KeyMax*urng();
                def_.keyR1 = ut::KeyMax*urng();
                def_.keyL2 = ut::KeyMax*urng();
                def_.keyR2 = ut::KeyMax*urng();
                
                return true;
            };
            
            Worm const& target(Worm const& worm) const {
                return worm;
            };
            
            template<typename Value>
            double ratio(data::Data<Value> const& data, state::State<Value>& state) const {
                return ut::beta()*ut::beta()*ut::beta()*ut::beta()
                    /((state.expansion()[def_.flavorL1].size() + 1)*(state.expansion()[def_.flavorR1].size() + 1)
                     *(state.expansion()[def_.flavorL2].size() + 1)*(state.expansion()[def_.flavorR2].size() + 1));
            };
        };
        
        
        template<typename Worm>
        struct QuadEraseStd : QuadErase {
            template<typename Value>
            QuadEraseStd(jsx::value const& jParams, data::Data<Value> const& data) : QuadErase(jParams, data) {};
            
            template<typename Value>
            bool propose(double const urn, data::Data<Value> const& data, state::State<Value>& state, ut::UniformRng& urng) {
                def_ = defs_[std::upper_bound(prob_.begin(), prob_.end(), urn*prob_.back()) - prob_.begin()];
                
                auto const& opsL1 = state.expansion()[def_.flavorL1];
                auto const& opsR1 = state.expansion()[def_.flavorR1];
                auto const& opsL2 = state.expansion()[def_.flavorL2];
                auto const& opsR2 = state.expansion()[def_.flavorR2];
                
                if(!(opsL1.size() && opsR1.size() && opsL2.size() && opsR2.size())) return false;
                if(!(opsL1 == opsL2 && opsL1.size()>1) or !(opsR1 == opsR2 && opsR1.size()>1)) return false;
                
                def_.keyL1 = opsL1[urng()*opsL1.size()];
                def_.keyR1 = opsR1[urng()*opsR1.size()];
                def_.keyL2 = opsL2[urng()*opsL2.size()];
                def_.keyR2 = opsR2[urng()*opsR2.size()];
                
                //Don't select repeated operators
                while (def_.keyL2 == def_.keyL1)
                    def_.keyL2 = opsL2[urng()*opsL2.size()];
                
                while (def_.keyR2 == def_.keyR1)
                    def_.keyR2 = opsR2[urng()*opsR2.size()];
                
                return true;
            };
            
            Worm const& target(Worm const& worm) const {
                return worm;
            };
            
            template<typename Value>
            double ratio(data::Data<Value> const& data, state::State<Value>& state) const {
                return (state.expansion()[def_.flavorL1].size()*state.expansion()[def_.flavorR1].size()
                       *state.expansion()[def_.flavorL2].size()*state.expansion()[def_.flavorR2].size())
                       /(ut::beta()*ut::beta()*ut::beta()*ut::beta());
            };
            
        };
        
    }
    
}

#endif
