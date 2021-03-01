#ifndef CTQMC_INCLUDE_UPDATES_EXPANSION_INSERTERASEQUAD_UPDATE_H
#define CTQMC_INCLUDE_UPDATES_EXPANSION_INSERTERASEQUAD_UPDATE_H


#include "../include/Generic.h"
#include "../../bath/algebra/Insert.h"
#include "../../bath/algebra/Erase.h"
#include "../../bath/algebra/MultiInsert.h"
#include "../../bath/algebra/MultiErase.h"
#include "../../bath/algebra/MultiInsert.h"
#include "../../bath/algebra/MultiErase.h"
#include "../../bath/algebra/MultiInsert.h"
#include "../../bath/algebra/MultiErase.h"


namespace upd {
    
    namespace expansion {
        
        
        struct QuadBase {
            QuadBase() = delete;
            template<typename Value>
            QuadBase(jsx::value const& jParams, data::Data<Value> const& data) {
                
                int bath1 = 0;
                for(auto const& block1 : data.hyb().blocks()) {
                    for(auto const& flavorL1 : block1.flavorsL())
                        for(auto const& flavorR1 : block1.flavorsR()){
                            
                            int bath2=0;
                            for(auto const& block2 : data.hyb().blocks()){
                                for(auto const& flavorL2 : block2.flavorsL())
                                    for(auto const& flavorR2 : block2.flavorsR()) {
                                        
                                        prob_.push_back((prob_.size() ? prob_.back() : .0) + 1.);
                                        defs_.push_back({flavorL1, flavorR1, flavorL2, flavorR2, bath1, bath2});
                                        
                                    }
                                bath2++;
                            }
                            
                        }
                    ++bath1;
                }
                
            }
            ~QuadBase() = default;
            
            jsx::value json() const {
                return jsx::null_t();
            };
            
        protected:
            struct Def {
                Def() = default;
                Def(int flavorL1, int flavorR1, int flavorL2, int flavorR2, int bath1, int bath2) :
                flavorL1(flavorL1), flavorR1(flavorR1), bath1(bath1),
                flavorL2(flavorL2), flavorR2(flavorR2), bath2(bath2) { same_bath = bath1 == bath2; };
                int flavorL1, flavorR1, bath1;
                int flavorL2, flavorR2, bath2;
                bool same_bath;
                ut::KeyType keyL1, keyR1;
                ut::KeyType keyL2, keyR2;
            };
            
            std::vector<double> prob_;
            std::vector<Def> defs_;
            Def def_;
        };
        
        
        struct QuadInsert : QuadBase {
            template<typename Value>
            QuadInsert(jsx::value const& jParams, data::Data<Value> const& data) : QuadBase(jParams, data) {};
            
            template<typename Mode, typename Value>
            bool impurity(data::Data<Value> const& data, state::State<Value>& state) const {
                
                if(!state.product().insert(def_.keyR1, def_.flavorR1)) return false;
                if(!state.product().insert(def_.keyL1, def_.flavorL1)) return false;
                if(!state.product().insert(def_.keyR2, def_.flavorR2)) return false;
                if(!state.product().insert(def_.keyL2, def_.flavorL2)) return false;
                
                state.dyn().insert(def_.keyR1, def_.flavorR1);
                state.dyn().insert(def_.keyL1, def_.flavorL1);
                state.dyn().insert(def_.keyR2, def_.flavorR2);
                state.dyn().insert(def_.keyL2, def_.flavorL2);
                
                return true;
            };
            
            template<typename Value>
            void bath(data::Data<Value> const& data, std::vector<bath::Bath<Value>>& baths) const {
                if (def_.same_bath){
                    bath::MultiInsert upd;
                    upd.VecInsert.push_back({bath::Insert{def_.keyL1, def_.flavorL1, def_.keyR1, def_.flavorR1}});
                    upd.VecInsert.push_back({bath::Insert{def_.keyL2, def_.flavorL2, def_.keyR2, def_.flavorR2}});
                    baths[def_.bath1].add(data.hyb(),upd);
                    
                }else{
                    baths[def_.bath1].add(data.hyb(), bath::Insert{def_.keyL1, def_.flavorL1, def_.keyR1, def_.flavorR1});
                    baths[def_.bath2].add(data.hyb(), bath::Insert{def_.keyL2, def_.flavorL2, def_.keyR2, def_.flavorR2});
                }
            };
            
            template<typename Value>
            void accept(data::Data<Value> const& data, state::State<Value>& state) const {
                
                state.expansion()[def_.flavorL1].insert(def_.keyL1);
                state.expansion()[def_.flavorR1].insert(def_.keyR1);
                state.expansion()[def_.flavorL2].insert(def_.keyL2);
                state.expansion()[def_.flavorR2].insert(def_.keyR2);
            };
            
            template<typename Value>
            void reject(data::Data<Value> const& data, state::State<Value>& state) {
            };
        };
        
        
        struct QuadErase : QuadBase {
            template<typename Value>
            QuadErase(jsx::value const& jParams, data::Data<Value> const& data) : QuadBase(jParams, data) {};
            
            template<typename Mode, typename Value>
            bool impurity(data::Data<Value> const& data, state::State<Value>& state) const {
                
                state.product().erase(def_.keyL2);
                state.product().erase(def_.keyR2);
                state.product().erase(def_.keyL1);
                state.product().erase(def_.keyR1);
                
                state.dyn().erase(def_.keyL2, def_.flavorL2);
                state.dyn().erase(def_.keyR2, def_.flavorR2);
                state.dyn().erase(def_.keyL1, def_.flavorL1);
                state.dyn().erase(def_.keyR1, def_.flavorR1);
                
                return true;
            };
            
            template<typename Value>
            void bath(data::Data<Value> const& data, std::vector<bath::Bath<Value>>& baths) const {
                
                if (def_.same_bath){
                    bath::MultiErase upd;
                    upd.VecErase.push_back({bath::Erase{def_.keyL1, def_.keyR1}});
                    upd.VecErase.push_back({bath::Erase{def_.keyL2, def_.keyR2}});
                    baths[def_.bath1].add(data.hyb(),upd);
                    
                }else{
                    baths[def_.bath1].add(data.hyb(), bath::Erase{def_.keyL1, def_.keyR1});
                    baths[def_.bath2].add(data.hyb(), bath::Erase{def_.keyL2, def_.keyR2});
                }
                
                
            };
            
            template<typename Value>
            void accept(data::Data<Value> const& data, state::State<Value>& state) const {
                
                state.expansion()[def_.flavorL1].erase(def_.keyL1);
                state.expansion()[def_.flavorR1].erase(def_.keyR1);
                state.expansion()[def_.flavorL2].erase(def_.keyL2);
                state.expansion()[def_.flavorR2].erase(def_.keyR2);
            };
            
            template<typename Value>
            void reject(data::Data<Value> const& data, state::State<Value>& state) {
            };
            
        };
        
        
    }
    
}

#endif
