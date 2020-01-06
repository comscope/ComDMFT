#ifndef UPDATES_IMPURITY_H
#define UPDATES_IMPURITY_H

#include <iostream>
#include <random>

#include "Definition.h"
#include "../Utilities.h"
#include "../Data.h"
#include "../State.h"

namespace imp {
    
    namespace itf {
        
        template<typename Def>
        struct Update {
            virtual ut::Flag surviving(Def const&, data::Data const&, state::State&) = 0;
            virtual ut::Zahl ratio(Def const&, data::Data const&, state::State&) = 0;
            virtual ut::Flag decide(ut::Zahl const&, data::Data const&, state::State&, itf::Batcher&) = 0;
            virtual void accept(Def const&, data::Data const&, state::State&) = 0;
            virtual void reject(Def const&, data::Data const&, state::State&) = 0;
            virtual ~Update() = default;
        };
        
    };
    
    template<typename Def, typename Alloc> struct Update;
    
    template<typename Alloc>
    struct Update<upd::InsertTwo, Alloc> : itf::Update<upd::InsertTwo> {
        ut::Flag surviving(upd::InsertTwo const& u, data::Data const& data, state::State& state) {
            if(!state.product().insert(u.keyR(), u.flavorR, u.ptrR())) return ut::Flag::Reject;
            if(!state.product().insert(u.keyL(), u.flavorL, u.ptrL())) return ut::Flag::Reject;
            
            densityMatrix_ = DensityMatrix<Alloc>(state.product(), data.eig());
            
            return densityMatrix_.surviving(data.eig());
        };
        
        ut::Zahl ratio(upd::InsertTwo const& u, data::Data const& data, state::State& state) {
            state.dyn().insert(u.keyR(), u.flavorR%2);
            state.dyn().insert(u.keyL(), u.flavorL%2);
            return state.dyn().ratio();
        };
        
        ut::Flag decide(ut::Zahl const& x, data::Data const& data, state::State& state, itf::Batcher& batcher) {
            return densityMatrix_.decide(state.densityMatrix().Z()*x, state.product(), batcher);
        };
        
        void accept(upd::InsertTwo const& u, data::Data const& data, state::State& state) {
            state.product().accept();
            get<Alloc>(state.densityMatrix()) = std::move(densityMatrix_);
            state.dyn().accept();
        };
        
        void reject(upd::InsertTwo const& u, data::Data const& data, state::State& state) {
            state.product().reject();
            densityMatrix_ = DensityMatrix<Alloc>();
            state.dyn().reject();
        };
    private:
        DensityMatrix<Alloc> densityMatrix_;
    };
    
    
    template<typename Alloc>
    struct Update<upd::EraseTwo, Alloc> : itf::Update<upd::EraseTwo> {
        ut::Flag surviving(upd::EraseTwo const& u, data::Data const& data, state::State& state) {
            state.product().erase(u.keyL());
            state.product().erase(u.keyR());
            
            densityMatrix_ = DensityMatrix<Alloc>(state.product(), data.eig());
            
            return densityMatrix_.surviving(data.eig());
        };
        
        ut::Zahl ratio(upd::EraseTwo const& u, data::Data const& data, state::State& state) {
            state.dyn().erase(u.keyL(), u.flavorL%2);
            state.dyn().erase(u.keyR(), u.flavorR%2);
            return state.dyn().ratio();
        };
        
        ut::Flag decide(ut::Zahl const& ratio, data::Data const& data, state::State& state, itf::Batcher& batcher) {
            return densityMatrix_.decide(state.densityMatrix().Z()*ratio, state.product(), batcher);
        };
        
        void accept(upd::EraseTwo const& u, data::Data const& data, state::State& state) {
            state.product().accept();
            get<Alloc>(state.densityMatrix()) = std::move(densityMatrix_);
            state.dyn().accept();
        };
        
        void reject(upd::EraseTwo const& u, data::Data const& data, state::State& state) {
            state.product().reject();
            densityMatrix_ = DensityMatrix<Alloc>();
            state.dyn().reject();
        };
    private:
        DensityMatrix<Alloc> densityMatrix_;
    };
}

#endif


