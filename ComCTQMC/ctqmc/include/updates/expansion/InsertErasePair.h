#ifndef CTQMC_INCLUDE_UPDATES_EXPANSION_INSERTERASEPAIR_UPDATE_H
#define CTQMC_INCLUDE_UPDATES_EXPANSION_INSERTERASEPAIR_UPDATE_H


#include "../include/Generic.h"
#include "../../bath/algebra/Insert.h"
#include "../../bath/algebra/Erase.h"


namespace upd {
    
    namespace expansion {
        
        
        struct Base {
            Base() = delete;
            template<typename Value>
            Base(jsx::value const& jParams, data::Data<Value> const& data) {
                int bath = 0;
                for(auto const& block : data.hyb().blocks()) {
                    for(auto const& flavorL : block.flavorsL())
                        for(auto const& flavorR : block.flavorsR()) {
                            prob_.push_back((prob_.size() ? prob_.back() : .0) + 1.);
                            defs_.push_back({flavorL, flavorR, bath});
                        }
                    ++bath;
                }
            }
            ~Base() = default;
            
            jsx::value json() const {
                return jsx::null_t();
            };
            
        protected:
            struct Def {
                Def() = default;
                Def(int flavorL, int flavorR, int bath) : flavorL(flavorL), flavorR(flavorR), bath(bath) {};
                int flavorL, flavorR, bath;
                ut::KeyType keyL, keyR;
            };
            
            std::vector<double> prob_;
            std::vector<Def> defs_;
            Def def_;
        };
        
        
        struct Insert : Base {
            template<typename Value>
            Insert(jsx::value const& jParams, data::Data<Value> const& data) : Base(jParams, data) {};
            
            template<typename Mode, typename Value>
            bool impurity(data::Data<Value> const& data, state::State<Value>& state) const {
                if(!state.product().insert(def_.keyR, def_.flavorR)) return false;
                if(!state.product().insert(def_.keyL, def_.flavorL)) return false;
                
                state.dyn().insert(def_.keyR, def_.flavorR);
                state.dyn().insert(def_.keyL, def_.flavorL);
                
                return true;
            };
            
            template<typename Value>
            void bath(data::Data<Value> const& data, std::vector<bath::Bath<Value>>& baths) const {
                baths[def_.bath].add(data.hyb(), bath::Insert{def_.keyL, def_.flavorL, def_.keyR, def_.flavorR});
            };
            
            template<typename Value>
            void accept(data::Data<Value> const& data, state::State<Value>& state) const {
                state.expansion()[def_.flavorL].insert(def_.keyL);
                state.expansion()[def_.flavorR].insert(def_.keyR);
            };
            
            template<typename Value>
            void reject(data::Data<Value> const& data, state::State<Value>& state) {
            };
        };
        
        
        struct Erase : Base {
            template<typename Value>
            Erase(jsx::value const& jParams, data::Data<Value> const& data) : Base(jParams, data) {};
            
            template<typename Mode, typename Value>
            bool impurity(data::Data<Value> const& data, state::State<Value>& state) const {
                state.product().erase(def_.keyL);
                state.product().erase(def_.keyR);
                
                state.dyn().erase(def_.keyL, def_.flavorL);
                state.dyn().erase(def_.keyR, def_.flavorR);
                
                return true;
            };
            
            template<typename Value>
            void bath(data::Data<Value> const& data, std::vector<bath::Bath<Value>>& baths) const {
                baths[def_.bath].add(data.hyb(), bath::Erase{def_.keyL, def_.keyR});
            };
            
            template<typename Value>
            void accept(data::Data<Value> const& data, state::State<Value>& state) const {
                state.expansion()[def_.flavorL].erase(def_.keyL);
                state.expansion()[def_.flavorR].erase(def_.keyR);
            };
            
            template<typename Value>
            void reject(data::Data<Value> const& data, state::State<Value>& state) {
            };
            
        };
        
        
    }
    
}

#endif
