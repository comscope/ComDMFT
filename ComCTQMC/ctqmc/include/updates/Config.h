#ifndef UPDATES_CONFIG_H
#define UPDATES_CONFIG_H

#include <stdexcept>
#include <iostream>
#include <vector>

#include "Definition.h"
#include "../Utilities.h"
#include "../Data.h"
#include "../State.h"

// Test low order !!!!! k = 0, 2 !!!!!!

namespace cfg {
    //----------------------------------------------------------------------------------------------------------------------------
    template<typename Def> struct Update;
    
    template<>
    struct Update<upd::InsertTwo> {
        int propose(upd::InsertTwo& u, data::Data const& data, state::State const& state, ut::UniformRng& urng) {
            auto const& opsL = state.config()[u.flavorL];
            auto const& opsR = state.config()[u.flavorR];
            
            ut::KeyType keyHigh;
            if(opsL.size() && opsR.size()) {
                std::size_t const index = urng()*(opsL.size() + opsR.size());
                ut::KeyType const keyLow = index < opsL.size() ? opsL[index].key() : opsR[index - opsL.size()].key();
                
                auto it = std::upper_bound(opsL.begin(), opsL.end(), Entry(keyLow));
                keyHigh = it != opsL.end() ? it->key() : opsL.begin()->key() + ut::KeyMax;
                
                it = std::upper_bound(opsR.begin(), opsR.end(), Entry(keyLow));
                keyHigh = std::min(it != opsR.end() ? it->key() : opsR.begin()->key() + ut::KeyMax, keyHigh);
                
                keyDiff_ = keyHigh - keyLow; if(keyHigh > ut::KeyMax) keyHigh -= ut::KeyMax;
            } else
                keyHigh = keyDiff_ = ut::KeyMax;
            
            
            entryL_ = Entry(ut::cyclic(keyHigh - urng()*keyDiff_));
            entryR_ = Entry(ut::cyclic(keyHigh - urng()*keyDiff_));
            u.keyL() = entryL_.key(); u.ptrL() = entryL_.ptr();
            u.keyR() = entryR_.key(); u.ptrR() = entryR_.ptr();

            return 1;
        };
        double ratio(upd::InsertTwo const& u, data::Data const& data, state::State const& state) {
            int const N = state.config()[u.flavorL].size() + state.config()[u.flavorR].size();
            double const timeDiff = (ut::beta()*keyDiff_)/ut::KeyMax;
            return timeDiff*timeDiff*(N ? N/(N + 2.) : 1.);
        };
        void accept(upd::InsertTwo const& u, data::Data const& data, state::State& state) {
            state.config()[u.flavorL].insert(std::move(entryL_));
            state.config()[u.flavorR].insert(std::move(entryR_));
        };
        void reject(upd::InsertTwo const& u, data::Data const& data, state::State& state) {
        };
        
    private:
        ut::KeyType keyDiff_;
        Entry entryL_;
        Entry entryR_;
    };
    
    template<>
    struct Update<upd::EraseTwo> {
        int propose(upd::EraseTwo& u, data::Data const& data, state::State const& state, ut::UniformRng& urng) {
            auto const& opsL = state.config()[u.flavorL];
            auto const& opsR = state.config()[u.flavorR];
            
            if(!(opsL.size() && opsR.size())) return 0;
            
            if(opsL.size() + opsR.size() > 2) {
                std::size_t index = urng()*(opsL.size() + opsR.size());
                ut::KeyType const keyLow = index < opsL.size() ? opsL[index].key() : opsR[index - opsL.size()].key();
                
                auto itL = std::upper_bound(opsL.begin(), opsL.end(), Entry(keyLow)); ut::KeyType shiftL = 0;
                if(itL != opsL.end()) u.keyL() = itL->key(); else { u.keyL() = (itL = opsL.begin())->key(); shiftL = ut::KeyMax;};
                
                auto itR = std::upper_bound(opsR.begin(), opsR.end(), Entry(keyLow)); ut::KeyType shiftR = 0;
                if(itR != opsR.end()) u.keyR() = itR->key(); else { u.keyR() = (itR = opsR.begin())->key(); shiftR = ut::KeyMax;};
                
                ut::KeyType keyHigh = ++itL != opsL.end() ? itL->key() + shiftL : opsL.begin()->key() + ut::KeyMax;
                keyHigh = std::min(++itR != opsR.end() ? itR->key() + shiftR : opsR.begin()->key() + ut::KeyMax, keyHigh);
                
                if(std::max(u.keyL() + shiftL, u.keyR() + shiftR) < keyHigh) keyDiff_ = keyHigh - keyLow; else return 0;
            } else {
                u.keyL() = opsL.begin()->key();
                u.keyR() = opsR.begin()->key();
                keyDiff_ = ut::KeyMax;
            }

            return 1;
        };
        double ratio(upd::EraseTwo const& u, data::Data const& data, state::State const& state) {
            int const N = state.config()[u.flavorL].size() + state.config()[u.flavorR].size();
            double const timeDiff = (ut::beta()*keyDiff_)/ut::KeyMax;
            return (N > 2 ? N/(N - 2.) : 1.)/(timeDiff*timeDiff);
        };
        void accept(upd::EraseTwo const& u, data::Data const& data, state::State& state) {
            state.config()[u.flavorL].erase(u.keyL());
            state.config()[u.flavorR].erase(u.keyR());
        };
        void reject(upd::EraseTwo const& u, data::Data const& data, state::State& state) {
        };
        
    private:
        ut::KeyType keyDiff_;
    };
    
}

#endif