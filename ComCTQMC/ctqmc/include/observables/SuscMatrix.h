#ifndef OBSERVABLES_SUSCMATRIX_H
#define OBSERVABLES_SUSCMATRIX_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

#include "Observable.h"
#include "../Utilities.h"
#include "../Data.h"
#include "../State.h"
#include "../../../include/measurements/Measurements.h"

namespace obs {
    
    template<typename Alloc, bool Direct, bool Bulla>
    struct SuscMatrix  : Observable {
        SuscMatrix() = delete;
        SuscMatrix(jsx::value const& jParams, data::Data const& data) :
        flavors_(data.ops().flavors()/2),
        nMat_(std::max(static_cast<int>(ut::beta()*jParams("susceptibility cutoff").real64()/(2*M_PI)), 1)),
        phase_(Phase::Calculate),
        urng_(std::mt19937(234), std::uniform_real_distribution<double>(.0, 1.)),
        accBulla_(flavors_, std::vector<std::vector<double>>(flavors_, std::vector<double>(nMat_, .0))),
        accDirect_(flavors_, std::vector<std::vector<double>>(flavors_, std::vector<double>(nMat_, .0))) {
        };
        SuscMatrix(SuscMatrix const&) = delete;
        SuscMatrix(SuscMatrix&&) = delete;
        SuscMatrix& operator=(SuscMatrix const&) = delete;
        SuscMatrix& operator=(SuscMatrix&&) = delete;
        ~SuscMatrix() = default;
        
        bool sample(ut::complex const sign, data::Data const& data, state::State& state, imp::itf::Batcher& batcher) {
            switch(phase_) {
                case Phase::Calculate:
                    calculate(data, imp::get<Alloc>(state.densityMatrix()), imp::get<Alloc>(state.product()), batcher);
                    
                    phase_ = Phase::Finalize; return false;
                    
                case Phase::Finalize:
                    finalize(sign, data, state);
                    
                    phase_ = Phase::Calculate;
            }
            
            return true;
        };
        
        void store(data::Data const& data, jsx::value& measurements, std::int64_t samples) {
            for(int f2 = 0; f2 < flavors_; ++f2)
                for(int f1 = 0; f1 < flavors_; ++f1) {
                    if(Bulla) {
                        measurements["susceptibility bulla"][std::to_string(f2) + "_" + std::to_string(f1)] << meas::fix(accBulla_[f2][f1], samples);
                        for(auto& x : accBulla_[f2][f1]) x = .0;
                    }
                    
                    if(Direct) {
                        measurements["susceptibility direct"][std::to_string(f2) + "_" + std::to_string(f1)] << meas::fix(accDirect_[f2][f1], samples);
                        for(auto& x : accDirect_[f2][f1]) x = .0;
                    }
                }
        };
        
        void finalize(data::Data const& data, jsx::value& measurements, std::int64_t samples) {
            store(data, measurements, samples);
        };
        
    private:
        
        enum class Phase { Calculate, Finalize };
        
        int const flavors_;
        std::size_t const nMat_;
        
        Phase phase_;
        ut::RandomNumberGenerator<std::mt19937, std::uniform_real_distribution<double> > urng_;  //should be removed in future ...
        double key2_, key1_;
        std::vector<ut::Zahl> valBTau2_, valBTau1_;
        std::vector<std::vector<ut::Zahl>> valBTau21_, valDTau21_;
        std::vector<std::vector<std::vector<double>>> accBulla_, accDirect_;
        
        
        void calculate(data::Data const& data, imp::DensityMatrix<Alloc> const& densityMatrix, imp::Product<Alloc>& product, imp::itf::Batcher& batcher) {
            auto const level = product.height();
            
            typename imp::Product<Alloc>::CAccessType node0 = product.first(), nodeA, nodeB;
            do { nodeA = product.insert(key1_ = urng_()*ut::KeyMax, &data.ide(), level + 1);} while(nodeA == product.last());
            do { nodeB = product.insert(key2_ = urng_()*ut::KeyMax, &data.ide(), level + 1);} while(nodeB == product.last());
            
            if(key2_ < key1_) std::swap(nodeA, nodeB);
            
            for(auto const sec0 : densityMatrix) {
                product.multiply(node0, level, sec0, batcher);
                auto const& mat0 = node0->op(level)->mat(sec0);
                
                int const secA = node0->op(level)->map(sec0).sector;
                product.multiply(nodeA, level, secA, batcher);
                auto const& matA = nodeA->op(level)->mat(secA);
                
                int const secB = nodeA->op(level)->map(secA).sector;
                product.multiply(nodeB, level, secB, batcher);
                auto const& matB = nodeB->op(level)->mat(secB);
                
                imp::mult(product.bufferC(), mat0, matB, batcher);
                
                int const sec2 = key2_ < key1_ ? secA : secB;
                int const sec1 = key2_ < key1_ ? secB : secA;
                
                auto const& mat2 = key2_ < key1_ ? matA : product.bufferC();
                auto const& mat1 = key2_ < key1_ ? product.bufferC() : matA;
                
                // [sec1,mat2,sec2]-key2-[sec2,mat1,sec1]-key1  if key2 > key1 and [sec2,mat1,sec1]-key1-[sec1,mat2,sec2]-key2 if key1 > key2.
                // We do not need to distinguish in the following since trace is cyclic.
                
                if(Bulla)
                    calculate_bulla(data, sec2, mat2, sec1, mat1, product.bufferA(), product.bufferB(), batcher);
                if(Direct)
                    calculate_direct(data, sec2, mat2, sec1, mat1, product.bufferA(), product.bufferB(), batcher);
            }
        };
        
        void calculate_bulla(data::Data const& data, int const sec2, imp::Matrix<Alloc> const& mat2, int const sec1, imp::Matrix<Alloc> const& mat1, imp::Matrix<Alloc>& bufferA, imp::Matrix<Alloc>& bufferB, imp::itf::Batcher& batcher) {
            auto const& bullaOcc = get<Alloc>(data.bullaOcc());
            
            valBTau1_.resize(flavors_, .0);
            valBTau2_.resize(flavors_, .0);
            valBTau21_.resize(flavors_, std::vector<ut::Zahl>(flavors_, .0));
            
            imp::mult(bufferA, mat2, mat1, batcher);
            for(int flavor1 = 0; flavor1 < flavors_; ++flavor1)
                if(bullaOcc.at(flavor1).map(sec1).sector != 0) {
                    auto const& bullaOcc1 = bullaOcc.at(flavor1).mat(sec1);
                    imp::traceAtB(nullptr, &valBTau1_[flavor1], bullaOcc1, bufferA, batcher);
                }
            
            for(int flavor2 = 0; flavor2 < flavors_; ++flavor2)
                if(bullaOcc.at(flavor2).map(sec2).sector != 0) {
                    auto const& bullaOcc2 = bullaOcc.at(flavor2).mat(sec2);
                    
                    imp::mult(bufferA, bullaOcc2, mat1, batcher);
                    imp::mult(bufferB, mat2, bufferA, batcher);
                    
                    imp::trace(nullptr, &valBTau2_[flavor2], bufferB, batcher);

                    for(int flavor1 = 0; flavor1 < flavors_; ++flavor1)
                        if(bullaOcc.at(flavor1).map(sec1).sector != 0) {
                            auto const& bullaOcc1 = bullaOcc.at(flavor1).mat(sec1);
                            imp::traceAtB(nullptr, &valBTau21_[flavor2][flavor1], bullaOcc1, bufferB, batcher);
                        }
                }
        };
        
        void calculate_direct(data::Data const& data, int const sec2, imp::Matrix<Alloc> const& mat2, int const sec1, imp::Matrix<Alloc> const& mat1, imp::Matrix<Alloc>& bufferA, imp::Matrix<Alloc>& bufferB, imp::itf::Batcher& batcher) {
            auto const& occ = get<Alloc>(data.occ());
            
            valDTau21_.resize(flavors_, std::vector<ut::Zahl>(flavors_, .0));
            for(int flavor2 = 0; flavor2 < flavors_; ++flavor2)
                if(occ.at(flavor2).map(sec2).sector != 0) {
                    auto const& occ2 = occ.at(flavor2).mat(sec2);
                    
                    imp::mult(bufferA, occ2, mat1, batcher);
                    imp::mult(bufferB, mat2, bufferA, batcher);
                    
                    for(int flavor1 = 0; flavor1 < flavors_; ++flavor1)
                        if(occ.at(flavor1).map(sec1).sector != 0) {
                            auto const& occ1 = occ.at(flavor1).mat(sec1);
                            imp::traceAtB(nullptr, &valDTau21_[flavor2][flavor1], occ1, bufferB, batcher);
                        }
                }
        };
        
        
        void finalize(ut::complex const sign, data::Data const& data, state::State& state) {
            state.product().reject();
            
            double const u2 = key2_/static_cast<double>(ut::KeyMax); ut::complex const fact2 = ut::complex(std::cos(2*M_PI*u2), std::sin(2*M_PI*u2));
            double const u1 = key1_/static_cast<double>(ut::KeyMax); ut::complex const fact1 = ut::complex(std::cos(2*M_PI*u1), std::sin(2*M_PI*u1));
            ut::complex exp2(1.), exp1(1.);
            
            std::vector<ut::complex> expTau1(nMat_), expTau2(nMat_);
            for(std::size_t n = 1; n < nMat_; ++n) {
                exp2 *= fact2; expTau2[n] = exp2;
                exp1 *= fact1; expTau1[n] = exp1;
            }
            
            if(Bulla)
                finalize_bulla(sign, expTau2, expTau1, state);
            if(Direct)
                finalize_direct(sign, expTau2, expTau1, state);
        };
    
        void finalize_bulla(ut::complex const sign, std::vector<ut::complex> const& expTau2, std::vector<ut::complex> const& expTau1, state::State& state) {
            std::vector<double> valTau1(flavors_), valTau2(flavors_);
            std::vector<std::vector<double>> valTau21(flavors_, std::vector<double>(flavors_));
            
            for(int f = 0; f < flavors_; ++f) valTau2[f] = (valBTau2_[f]/state.densityMatrix().Z()).to_double();
            for(int f = 0; f < flavors_; ++f) valTau1[f] = -(valBTau1_[f]/state.densityMatrix().Z()).to_double(); //gets minus sign because bulla^\dagger = -bulla !!
            for(int f2 = 0; f2 < flavors_; ++f2)
                for(int f1 = 0; f1 < flavors_; ++f1)
                    valTau21[f2][f1] = -(valBTau21_[f2][f1]/state.densityMatrix().Z()).to_double();               //gets minus sign because bulla^\dagger = -bulla !!
            
            valBTau2_.clear(); valBTau1_.clear(); valBTau21_.clear();
            
            std::vector<double> zeroTau2(flavors_, .0), zeroTau1(flavors_, .0);
            
            for(int f = 0; f < flavors_; ++f) {
                auto& entry2 = zeroTau2[f];
                for(auto const& c : state.config()[2*f]) {
                    double const time = ut::beta()*std::abs((c.key() - key2_)/static_cast<double>(ut::KeyMax));
                    entry2 += time*(ut::beta() - time);
                }
                for(auto const& c : state.config()[2*f + 1]) {
                    double const time = ut::beta()*std::abs((c.key() - key2_)/static_cast<double>(ut::KeyMax));
                    entry2 -= time*(ut::beta() - time);
                }
                
                auto& entry1 = zeroTau1[f];
                for(auto const& c : state.config()[2*f]) {
                    double const time = ut::beta()*std::abs((c.key() - key1_)/static_cast<double>(ut::KeyMax));
                    entry1 += time*(ut::beta() - time);
                }
                for(auto const& c : state.config()[2*f + 1]) {
                    double const time = ut::beta()*std::abs((c.key() - key1_)/static_cast<double>(ut::KeyMax));
                    entry1 -= time*(ut::beta() - time);
                }
            }

            double const time = ut::beta()*std::abs((key1_ - key2_)/static_cast<double>(ut::KeyMax));
            for(int f2 = 0; f2 < flavors_; ++f2)
                for(int f1 = 0; f1 < flavors_; ++f1) {
                    auto& entry = accBulla_[f2][f1];
                    
                    entry[0] += sign.real()*(
                                             +(time*(ut::beta() - time)*valTau21[f2][f1])*ut::beta()
                                             -(valTau2[f2]*zeroTau2[f1] + zeroTau2[f2]*valTau2[f1] +
                                               valTau1[f2]*zeroTau1[f1] + zeroTau1[f2]*valTau1[f1])/2.
                                             );
                    
                    auto const& expFlavor2 = state.exponentials().at(f2);
                    auto const& expFlavor1 = state.exponentials().at(f1);
                    
                    // TODO: sample negative matsubara frequencies as well if hybridisation is complex
                    for(std::size_t n = 1; n < nMat_; ++n)
                        entry[n] += (sign*(
                                          +(expTau2[n]*std::conj(expTau1[n])*valTau21[f2][f1])*ut::beta()
                                          -((expTau2[n]*valTau2[f2] + expTau1[n]*valTau1[f2])*std::conj(expFlavor1[n]) +
                                             expFlavor2[n]*std::conj(expTau2[n]*valTau2[f1] + expTau1[n]*valTau1[f1]))/2.
                                          )
                                     ).real();
                }
            
        };
        
        void finalize_direct(ut::complex const sign, std::vector<ut::complex> const& expTau2, std::vector<ut::complex> const& expTau1, state::State& state) {
            std::vector<std::vector<double>> valTau21(flavors_, std::vector<double>(flavors_));
            
            for(int f2 = 0; f2 < flavors_; ++f2)
                for(int f1 = 0; f1 < flavors_; ++f1)
                    valTau21[f2][f1] = (valDTau21_[f2][f1]/state.densityMatrix().Z()).to_double();
            
            valDTau21_.clear();
            
            for(int f2 = 0; f2 < flavors_; ++f2)
                for(int f1 = 0; f1 < flavors_; ++f1) {
                    auto& entry = accDirect_[f2][f1];
                    
                    entry[0] += sign.real()*valTau21[f2][f1]*ut::beta();
                    
                    // TODO: sample negative matsubara frequencies as well if hybridisation is complex
                    for(std::size_t n = 1; n < nMat_; ++n)
                        entry[n] += (sign*expTau2[n]*std::conj(expTau1[n])*valTau21[f2][f1]).real()*ut::beta();
                }
        };
        
    };
}

#endif
