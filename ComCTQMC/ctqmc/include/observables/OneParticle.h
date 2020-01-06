#ifndef OBSERVABLES_ONEPARTICLE_H
#define OBSERVABLES_ONEPARTICLE_H

#include <cmath>
#include <stdexcept>
#include <iostream>

#include "Observable.h"
#include "../Utilities.h"
#include "../Data.h"
#include "../State.h"
#include "../../../include/measurements/Measurements.h"

namespace obs {
    
    struct Green {
        static std::string name() { return "green";};
        static double value(ut::complex sign, double value, double* ptrL, double* ptrR) { return sign.real()*value;};
        static ut::complex value(ut::complex sign, ut::complex value, double* ptrL, double* ptrR) { return sign*value;};
    };
    struct BullaL {
        static std::string name() { return "bullaL";};
        static double value(ut::complex sign, double value, double* ptrL, double* ptrR) { return sign.real()*value**ptrL;};
        static ut::complex value(ut::complex sign, ut::complex value, double* ptrL, double* ptrR) { return sign*value**ptrL;};
    };
    struct BullaR {
        static std::string name() { return "bullaR";};
        static double value(ut::complex sign, double value, double* ptrL, double* ptrR) { return sign.real()*value**ptrR;};
        static ut::complex value(ut::complex sign, ut::complex value, double* ptrL, double* ptrR) { return sign*value**ptrR;};
    };
    

    template<template<typename> class Meas, typename Select, typename HybVal>
    struct OneParticle : Observable {
        OneParticle() = delete;
        OneParticle(jsx::value const& jParams, data::Data const& data) :
        flavors_(2*jParams("hybridisation")("matrix").size()),
        matrix_(flavors_*flavors_, nullptr)
        {
            auto const& jMatrix = jParams("hybridisation")("matrix");
            for(std::size_t i = 0; i < jMatrix.size(); ++i)
                for(std::size_t j = 0; j < jMatrix.size(); ++j)
                    if(jMatrix(i)(j).string() != "") {
                        auto const& entry = jMatrix(i)(j).string();
                        
                        if(!data_.count(entry))
                            data_.emplace(entry, std::make_pair(0, Meas<HybVal>(jParams)));
                        
                        data_.at(entry).first += 1;
                        matrix_[2*j  + flavors_*(2*i + 1)] = &data_.at(entry).second;  //hyb isch symmetrisch wil reell, c.f. Hyb.h Muscht du aber nochmal gucken wegen allgemeinem fall ....
                    }
        };
        OneParticle(OneParticle const&) = delete;
        OneParticle(OneParticle&&) = delete;
        OneParticle& operator=(OneParticle const&) = delete;
        OneParticle& operator=(OneParticle&&) = delete;
        ~OneParticle() = default;
        
        bool sample(ut::complex const sign, data::Data const& data, state::State& state, imp::itf::Batcher& batcher) {
            for(int b = 0; b < data.hyb().blocks().size(); ++b) {
                auto const& bath = bath::get<HybVal>(state.bath(b));
                auto data = bath.B().data();
                for(auto const& opL : bath.opsL())
                    for(auto const& opR : bath.opsR())
                        matrix_[opR.flavor() + flavors_*opL.flavor()]->add(opR.key() - opL.key(), Select::value(sign, *data++, opL.ptr(), opR.ptr()));
            }

            return true;
        };
        
        void store(data::Data const& data, jsx::value& measurements, std::int64_t samples) {
            for(auto& entry : data_)
                entry.second.second.store(measurements[Select::name()][entry.first], entry.second.first*samples);
        };
        
        void finalize(data::Data const& data, jsx::value& measurements, std::int64_t samples) {
            store(data, measurements, samples);
        };

    private:
        int const flavors_;
        
        std::map<std::string, std::pair<std::int64_t, Meas<HybVal>>> data_;
        std::vector<Meas<HybVal>*> matrix_;
    };
    
    
    template<typename HybVal>
    struct Vec4 {
        Vec4() = default;
        Vec4(HybVal x) : at{x, x, x, x} {};
        Vec4(Vec4 const&) = default;
        Vec4(Vec4&&) = default;
        Vec4& operator=(Vec4 const&) = default;
        Vec4& operator=(Vec4&&) = default;
        ~Vec4() = default;
        
        HybVal at[4];
        
        void operator+=(Vec4 const& other) {
            at[0] += other.at[0]; at[1] += other.at[1]; at[2] += other.at[2]; at[3] += other.at[3];
        };
        void operator*=(Vec4<double> const& other) {
            at[0] *= other.at[0]; at[1] *= other.at[1]; at[2] *= other.at[2]; at[3] *= other.at[3];
        };
        HybVal reduce() const {
            return at[0] + at[1] + at[2] + at[3];
        };
    };
    
    template<typename HybVal>
    inline Vec4<HybVal> operator+(Vec4<HybVal> const& lhs, Vec4<HybVal> const& rhs) {
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
    
    
    template<typename HybVal>
    struct Legendre {
        Legendre() = delete;
        Legendre(jsx::value const& jParams) :
        nPol_(jParams("green legendre cutoff").int64()),
        recursion_(nPol_),
        pos_(0),
        coeff_(nPol_, Vec4<HybVal>(.0))
        {
            for(std::size_t n = 2; n < nPol_; ++n) recursion_[n] = std::make_pair(Vec4<double>((2.*n - 1.)/n), Vec4<double>(-(n - 1.)/n));
        };
        Legendre(Legendre const&) = delete;
        Legendre(Legendre&&) = default;
        Legendre& operator=(Legendre const&) = delete;
        Legendre& operator=(Legendre&&) = delete;
        ~Legendre() = default;
        
        void add(ut::KeyType key, HybVal value) {
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
            
            std::vector<HybVal> coeff(nPol_);
            for(std::size_t n = 0; n < coeff_.size(); ++n)
                coeff[n] = -(2.*n + 1)/ut::beta()*coeff_[n].reduce();   //missing -1/beta factor

            measurements << meas::fix(coeff, samples);
            
            std::fill(coeff_.begin(), coeff_.end(), Vec4<HybVal>(.0));
        };
        
    private:
        std::size_t const nPol_;
        std::vector<std::pair<Vec4<double>, Vec4<double>>> recursion_;
        
        std::size_t pos_;
        Vec4<double> x_; Vec4<HybVal> val_;
        std::vector<Vec4<HybVal>> coeff_;

        void vec_add() {
            Vec4<HybVal> p0 = val_; coeff_[0] += p0;
            Vec4<HybVal> p1 = x_*val_; coeff_[1] += p1;
            
            for(std::size_t n = 2; n < nPol_; ++n) {
                Vec4<HybVal> p = (recursion_[n].first*x_)*p1 + recursion_[n].second*p0;
                coeff_[n] += p;
                p0 = p1; p1 = p;
            }
            
            pos_ = 0;
        };
    };
    
    template<typename HybVal>
    struct BinMoments {
        BinMoments() = delete;
        BinMoments(jsx::value const& jParams) :
        nMatG_(std::max(static_cast<int>(ut::beta()*jParams("green matsubara cutoff").real64()/(2*M_PI)), 1)),
        nItG_(4*(2*nMatG_ + 1)),
        DeltaInv_(nItG_/ut::beta()),
        data_(4*nItG_, .0) {
        };
        BinMoments(BinMoments const&) = delete;
        BinMoments(BinMoments&&) = default;
        BinMoments& operator=(BinMoments const&) = delete;
        BinMoments& operator=(BinMoments&&) = delete;
        ~BinMoments() = default;
        
        void add(ut::KeyType key, HybVal value) {
            if(key < 0) {
                key += ut::KeyMax;
                value *= -1.;
            }
            
            double const time = key*ut::beta()/ut::KeyMax;
            int const index = static_cast<int>(DeltaInv_*time);
            double const Dtime = time - static_cast<double>(index + .5)/DeltaInv_;
            
            HybVal* data = data_.data() + 4*index;
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
        
        std::vector<HybVal> data_;
    };
}

#endif