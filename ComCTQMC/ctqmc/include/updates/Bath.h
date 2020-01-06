#ifndef UPDATES_BATH_H
#define UPDATES_BATH_H

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include "Definition.h"
#include "../Utilities.h"
#include "../Data.h"
#include "../State.h"

//Vielleicht insert und erase noch verändern wie in CTBI hürütütütüt pluschter pluschter

namespace bath {
    
    namespace itf {
        
        template<typename Def>
        struct Update {
            virtual double ratio(Def const&, data::Data const&, state::State const&) = 0;
            virtual void accept(Def const&, data::Data const&, state::State&) = 0;
            virtual void reject(Def const&, data::Data const&, state::State&) = 0;
            virtual ~Update() = default;
        };
        
    };

    template<typename Def, typename HybVal> struct Update;
    
    template<typename HybVal>
    struct Update<upd::InsertTwo, HybVal> : itf::Update<upd::InsertTwo> {
        double ratio(upd::InsertTwo const& u, data::Data const& data, state::State const& state) {
            auto const& bath = get<HybVal>(state.bath(u.bath));
            auto const& hyb = get<HybVal>(data.hyb());
            
            int const N = bath.opsL_.size();
            val_ = hyb(u.flavorL, u.flavorR, u.keyL() - u.keyR());
            
            if(N) {
                Bv_.resize(N); vec_.resize(N);
                for(int n = 0; n < N; ++n) vec_[n] = hyb(bath.opsL_[n].flavor(), u.flavorR, bath.opsL_[n].key() - u.keyR());
                
                char const no = 'n';
                int const inc = 1;
                HybVal const zero = .0;
                HybVal const one = 1.;
                gemv(&no, &N, &N, &one, bath.B_.data(), &N, vec_.data(), &inc, &zero, Bv_.data(), &inc);
                
                for(int n = 0; n < N; ++n) vec_[n] = hyb(u.flavorL, bath.opsR_[n].flavor(), u.keyL() - bath.opsR_[n].key());
                val_ -= xTy(&N, vec_.data(), &inc, Bv_.data(), &inc);
            }
            
            return std::abs(val_);
        };

        void accept(upd::InsertTwo const& u, data::Data const& data, state::State& state) {
            auto& bath = get<HybVal>(state.bath(u.bath));
            
            int const N = bath.opsL_.size();
            int const newN = N + 1;
            
            bath.posL_[u.keyL()] = bath.opsL_.size(); bath.opsL_.push_back(Operator(u.keyL(), u.flavorL, u.ptrL()));
            bath.posR_[u.keyR()] = bath.opsR_.size(); bath.opsR_.push_back(Operator(u.keyR(), u.flavorR, u.ptrR()));
            
            Matrix<HybVal> temp(newN);
            
            HybVal fact = 1./val_;
            temp.at(N, N) = fact;
            
            if(N) {
                std::vector<HybVal> hBTilde(N);
                
                char const yes = 't';
                int const inc = 1;
                HybVal const zero = .0;
                HybVal const one = 1.;
                gemv(&yes, &N, &N, &fact, bath.B_.data(), &N, vec_.data(), &inc, &zero, hBTilde.data(), &inc);
                axyT(&N, &N, &one, Bv_.data(), &inc, hBTilde.data(), &inc, bath.B_.data(), &N);
                
                for(int n = 0; n < N; ++n)
                    copy(&N, bath.B_.data(0, n), &inc, temp.data(0, n), &inc);
                
                fact = -1./val_;
                scal(&N, &fact, Bv_.data(), &inc);
                copy(&N, Bv_.data(), &inc, temp.data(0, N), &inc);
                
                HybVal const minus = -1.;
                scal(&N, &minus, hBTilde.data(), &inc);
                copy(&N, hBTilde.data(), &inc, temp.data(N, 0), &newN);
            }
            
            bath.B_ = std::move(temp);
            
            bath.det_ *= val_;
        };
        
        void reject(upd::InsertTwo const& u, data::Data const& data, state::State& state) {
        };
    private:
        std::vector<HybVal> Bv_;
        std::vector<HybVal> vec_;
        
        HybVal val_;
    };

    template<typename HybVal>
    struct Update<upd::EraseTwo, HybVal> : itf::Update<upd::EraseTwo> {
        double ratio(upd::EraseTwo const& u, data::Data const& data, state::State const& state) {
            auto const& bath = get<HybVal>(state.bath(u.bath));
            
            itL_ = bath.posL_.find(u.keyL()); itR_ = bath.posR_.find(u.keyR());
            
            return std::abs(val_ = bath.B_.at(itR_->second, itL_->second));
        };
        
        void accept(upd::EraseTwo const& u, data::Data const& data, state::State& state) {
            auto& bath = get<HybVal>(state.bath(u.bath));
            
            int const N = bath.opsL_.size(); int const newN = N - 1;
            int const posL = itL_->second; int const posR = itR_->second;
            
            Matrix<HybVal> temp(newN);
            
            if(newN) {
                int const inc = 1;
        
                if(posL != newN) {
                    swap(&N, bath.B_.data(0, newN), &inc, bath.B_.data(0, posL), &inc); bath.swapSign_ *= -1;
                    bath.posL_[bath.opsL_.back().key()] = posL; bath.opsL_[posL] = bath.opsL_.back();
                }
                if(posR != newN) {
                    swap(&N, bath.B_.data(newN, 0), &N, bath.B_.data(posR, 0), &N); bath.swapSign_ *= -1;
                    bath.posR_[bath.opsR_.back().key()] = posR; bath.opsR_[posR] = bath.opsR_.back();
                }
                
                for(int n = 0; n < newN; ++n)
                    copy(&newN, bath.B_.data(0, n), &inc, temp.data(0, n), &inc);
                
                HybVal const fact = -1./val_;
                axyT(&newN, &newN, &fact, bath.B_.data(0, newN), &inc, bath.B_.data(newN, 0), &N, temp.data(), &newN);
            }
            
            bath.B_ = std::move(temp);
            
            bath.posL_.erase(itL_); bath.opsL_.pop_back();
            bath.posR_.erase(itR_); bath.opsR_.pop_back();
            
            bath.det_ *= val_;
        };
        
        void reject(upd::EraseTwo const& u, data::Data const& data, state::State& state) {
        };
    private:
        HybVal val_;
        
        std::map<ut::KeyType, int>::const_iterator itL_;
        std::map<ut::KeyType, int>::const_iterator itR_;
    };
}

#endif
