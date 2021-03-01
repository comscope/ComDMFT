#ifndef CTQMC_INCLUDE_OBSERVABLES_PARTITION_ONEPARTICLE_H
#define CTQMC_INCLUDE_OBSERVABLES_PARTITION_ONEPARTICLE_H

#include <cmath>
#include <stdexcept>
#include <iostream>

#include "../Observable.h"
#include "../../Utilities.h"
#include "../../Data.h"
#include "../../State.h"
#include "../../../../include/measurements/Measurements.h"

namespace obs {
    
    namespace partition {
        
        struct Green {
            static std::string name() { return "green";};
            template<typename Value> static Value get(Value value, Value bullaL, Value bullaR) { return value;};
        };
        struct BullaL {
            static std::string name() { return "bullaL";};
            template<typename Value> static Value get(Value value, Value bullaL, Value bullaR) { return value*bullaL;};
        };
        struct BullaR {
            static std::string name() { return "bullaR";};
            template<typename Value> static Value get(Value value, Value bullaL, Value bullaR) { return value*bullaR;};
        };
        
        
        template<template<typename> class Meas, typename Select, typename Value>
        struct OneParticle : obs::itf::Observable<Value> {
            OneParticle() = delete;
            OneParticle(std::int64_t store, jsx::value const& jParams, data::Data<Value> const& data) :
            flavors_(2*jParams("hybridisation")("matrix").size()),
            store_(store), samples_(0),
            matrix_(flavors_*flavors_, nullptr)
            {
                auto const& jMatrix = jParams("hybridisation")("matrix");
                for(std::size_t i = 0; i < jMatrix.size(); ++i)
                    for(std::size_t j = 0; j < jMatrix.size(); ++j)
                        if(jMatrix(i)(j).string() != "") {
                            auto const& entry = jMatrix(i)(j).string();
                            
                            if(!data_.count(entry))
                                data_.emplace(entry, std::make_pair(0, Meas<Value>(jParams)));
                            
                            // conceptually, the indexing should be (j, i) since \Delta_ij gives \G_ji, but for evalsim, it is handier
                            // to have the Green function represented by the "hybridisation" : "matrix" entry in the params.json
                            data_.at(entry).first += 1;
                            matrix_[2*i  + flavors_*(2*j + 1)] = &data_.at(entry).second;
                        }
            };
            OneParticle(OneParticle const&) = delete;
            OneParticle(OneParticle&&) = delete;
            OneParticle& operator=(OneParticle const&) = delete;
            OneParticle& operator=(OneParticle&&) = delete;
            ~OneParticle() = default;
            
            bool sample(Value const sign, data::Data<Value> const& data, state::State<Value>& state, jsx::value& measurements, imp::itf::Batcher<Value>& batcher) {
                for(auto const& bath : state.baths()) {
                    auto data = bath.B().data();
                    for(auto const& opL : bath.opsL())
                        for(auto const& opR : bath.opsR())
                            matrix_[opR.flavor() + flavors_*opL.flavor()]->add(opR.key() - opL.key(), sign*Select::get(*data++ , opL.bulla(), opR.bulla()));
                }
                
                ++samples_; if(samples_%store_ == 0) store(data, measurements);
                
                return true;
            };
            
            void finalize(data::Data<Value> const& data, jsx::value& measurements) {
                store(data, measurements);
            };
            
        private:
            int const flavors_;
            
            std::int64_t const store_;
            std::int64_t samples_;
            
            std::map<std::string, std::pair<std::int64_t, Meas<Value>>> data_;
            std::vector<Meas<Value>*> matrix_;
            
            
            void store(data::Data<Value> const& data, jsx::value& measurements) {
                for(auto& entry : data_)
                    entry.second.second.store(measurements[Select::name()][entry.first], entry.second.first*samples_);
                
                samples_ = 0;
            };
        };
        
        
        template<typename Value>
        struct Vec4 {
            Vec4() = default;
            Vec4(Value x) : at{x, x, x, x} {};
            Vec4(Vec4 const&) = default;
            Vec4(Vec4&&) = default;
            Vec4& operator=(Vec4 const&) = default;
            Vec4& operator=(Vec4&&) = default;
            ~Vec4() = default;
            
            Value at[4];
            
            void operator+=(Vec4 const& other) {
                at[0] += other.at[0]; at[1] += other.at[1]; at[2] += other.at[2]; at[3] += other.at[3];
            };
            void operator*=(Vec4<double> const& other) {
                at[0] *= other.at[0]; at[1] *= other.at[1]; at[2] *= other.at[2]; at[3] *= other.at[3];
            };
            Value reduce() const {
                return at[0] + at[1] + at[2] + at[3];
            };
        };
        
        template<typename Value>
        inline Vec4<Value> operator+(Vec4<Value> const& lhs, Vec4<Value> const& rhs) {
            auto temp = lhs; temp += rhs; return temp;
        };
        inline Vec4<double> operator*(Vec4<double> const& lhs, Vec4<double> const& rhs) {
            auto temp = lhs; temp *= rhs; return temp;
        };
        inline Vec4<ut::complex> operator*(Vec4<double> const& lhs, Vec4<ut::complex> const& rhs) {
            auto temp = rhs; temp *= lhs; return temp;
        };
        inline Vec4<ut::complex> operator*(Vec4<ut::complex> const& lhs, Vec4<double> const& rhs) {
            auto temp = lhs; temp *= rhs; return temp;
        };
        
        
        template<typename Value>
        struct Legendre {
            Legendre() = delete;
            Legendre(jsx::value const& jParams) :
            nPol_(jParams(cfg::partition::Worm::name())("green legendre cutoff").int64()),
            recursion_(nPol_),
            pos_(0),
            coeff_(nPol_, Vec4<Value>(.0))
            {
                for(std::size_t n = 2; n < nPol_; ++n) recursion_[n] = std::make_pair(Vec4<double>((2.*n - 1.)/n), Vec4<double>(-(n - 1.)/n));
            };
            Legendre(Legendre const&) = delete;
            Legendre(Legendre&&) = default;
            Legendre& operator=(Legendre const&) = delete;
            Legendre& operator=(Legendre&&) = delete;
            ~Legendre() = default;
            
            void add(ut::KeyType key, Value value) {
                if(key < 0) {
                    key += ut::KeyMax;
                    value *= -1.;
                }
                x_.at[pos_] = (key - ut::KeyMax/2)*(2./ut::KeyMax); val_.at[pos_] = value;
                
                if(++pos_ == 4) vec_add();
            };
            
            void store(jsx::value& measurements, std::int64_t samples) {
                if(pos_) {
                    for(; pos_ < 4; ++pos_) { x_.at[pos_] = .0; val_.at[pos_] = .0;};
                    vec_add();
                };
                
                std::vector<Value> coeff(nPol_);
                for(std::size_t n = 0; n < coeff_.size(); ++n)
                    coeff[n] = -(2.*n + 1)/ut::beta()*coeff_[n].reduce();   //missing -1/beta factor
                
                measurements << meas::fix(coeff, samples);
                
                std::fill(coeff_.begin(), coeff_.end(), Vec4<Value>(.0));
            };
            
        private:
            std::size_t const nPol_;
            std::vector<std::pair<Vec4<double>, Vec4<double>>> recursion_;
            
            std::size_t pos_;
            Vec4<double> x_; Vec4<Value> val_;
            std::vector<Vec4<Value>> coeff_;
            
            void vec_add() {
                Vec4<Value> p0 = val_; coeff_[0] += p0;
                Vec4<Value> p1 = x_*val_; coeff_[1] += p1;
                
                for(std::size_t n = 2; n < nPol_; ++n) {
                    Vec4<Value> p = (recursion_[n].first*x_)*p1 + recursion_[n].second*p0;
                    coeff_[n] += p;
                    p0 = p1; p1 = p;
                }
                
                pos_ = 0;
            };
        };
        
        template<typename Value>
        struct BinMoments {
            BinMoments() = delete;
            BinMoments(jsx::value const& jParams) :
            nMatG_(std::max(static_cast<int>(ut::beta()*jParams(cfg::partition::Worm::name())("green matsubara cutoff").real64()/(2*M_PI)), 1)),
            nItG_(4*(2*nMatG_ + 1)),
            DeltaInv_(nItG_/ut::beta()),
            data_(4*nItG_, .0) {
            };
            BinMoments(BinMoments const&) = delete;
            BinMoments(BinMoments&&) = default;
            BinMoments& operator=(BinMoments const&) = delete;
            BinMoments& operator=(BinMoments&&) = delete;
            ~BinMoments() = default;
            
            void add(ut::KeyType key, Value value) {
                if(key < 0) {
                    key += ut::KeyMax;
                    value *= -1.;
                }
                
                double const time = key*ut::beta()/ut::KeyMax;
                int const index = static_cast<int>(DeltaInv_*time);
                double const Dtime = time - static_cast<double>(index + .5)/DeltaInv_;
                
                Value* data = data_.data() + 4*index;
                *data++ += value;
                value *= Dtime;
                *data++ += value;
                value *= Dtime;
                *data++ += value;
                value *= Dtime;
                *data += value;
            };
            
            void store(jsx::value& measurements, std::int64_t samples) {
                measurements << meas::fix(data_, samples);
                
                std::fill(data_.begin(), data_.end(), .0);
            };
            
        private:
            int const nMatG_;
            int const nItG_;
            double const DeltaInv_;
            
            std::vector<Value> data_;
        };
        
    }
}

#endif