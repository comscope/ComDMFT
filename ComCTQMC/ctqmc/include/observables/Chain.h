#ifndef OBSERVABLES_CHAIN_H
#define OBSERVABLES_CHAIN_H

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include "Observable.h"
#include "../Utilities.h"
#include "../Data.h"
#include "../State.h"
#include "../../../include/measurements/Measurements.h"

namespace obs {
    
    template<typename Alloc, typename HybVal, bool Density, bool Bulla>
    struct Chain : Observable {
        Chain() = delete;
        Chain(jsx::value const& jParams, data::Data const& data) :
        phase_(Phase::Prepare),
        samples_(0) {
            init(imp::get<Alloc>(data.eig()), HybVal());
        };
        Chain(Chain const&) = delete;
        Chain(Chain&&) = delete;
        Chain operator=(Chain const&) = delete;
        Chain& operator=(Chain&&) = delete;
        ~Chain() = default;
        
        bool sample(ut::complex const sign, data::Data const& data, state::State& state, imp::itf::Batcher& batcher) { // scheisse das ...
            auto const& densityMatrix = imp::get<Alloc>(state.densityMatrix());
            auto& product = imp::get<Alloc>(state.product());
            
            if(product.size()) {
                switch(phase_) {
                    case Phase::Prepare:
                        for(auto it = product.first(); it != product.last(); it = it.next(0))
                            chain_.emplace_back(it.key(), it->flavor, it->prop(), it->ptr);
                        
                        it_ = densityMatrix.begin();
                        
                        phase_ = Phase::Calculate;
                        
                    case Phase::Calculate:
                        if(it_ != densityMatrix.end()) {
                            evaluate_static(sign, *it_, data, state, batcher);
                            ++it_; return false;
                        }
                        
                        if(Bulla) {
                            for(auto it = chain_.begin() + 1; it != chain_.end(); ++it) *it->ptr = (it->acc/densityMatrix.Z()).to_double();
                            if(data.dyn() != nullptr) evaluate_dynamic(data, state);
                        }
                        
                        chain_.clear();
                        
                        phase_ = Phase::Prepare;
                }
            } else if(Density) {
                auto const& ide = imp::get<Alloc>(data.ide());
                
                for(auto const sec : densityMatrix) {
                    imp::density_matrix(product.bufferA(), ide.mat(sec), product.first()->prop()->at(sec), ide.mat(sec), batcher);
                    accDensityMatrix(sign, sec, 1./(ut::beta()*densityMatrix.Z()), product.bufferA(), batcher, HybVal());
                }
            }
            
            ++samples_; return true;
        };
        
        void store(data::Data const&, jsx::value&, std::int64_t) {};
        
        void finalize(data::Data const& data, jsx::value& measurements, std::int64_t) {
            if(Density) {
                if(!measurements["density matrix"].is<jsx::array_t>())
                    measurements["density matrix"] = jsx::array_t(data.eig().sectorNumber());
                
                for(int sec = 0; sec < data.eig().sectorNumber(); ++sec)
                    store(imp::get<Alloc>(data.eig()), measurements["density matrix"], sec, HybVal());
            }
        };
        
    private:
        
        enum class Phase { Prepare, Calculate };
        
        struct Entry {
            Entry() = delete;
            Entry(ut::KeyType key, int flavor, imp::Propagator<Alloc> const* prop, double* ptr) :
            key(key), flavor(flavor), prop(prop), acc(.0), ptr(ptr) {
            };
            Entry(Entry const&) = delete;
            Entry(Entry&&) = default;
            Entry& operator=(Entry const&) = delete;
            Entry& operator=(Entry&&) = default;
            ~Entry() = default;
            
            ut::KeyType const key;
            int const flavor;
            imp::Propagator<Alloc> const* const prop;
            ut::Zahl acc;
            double* const ptr;
            
            std::unique_ptr<imp::Matrix<Alloc>> rightDensity, rightBulla;
        };

    
        Phase phase_;
        std::int64_t samples_;
        
        std::vector<Entry> chain_;
        typename imp::DensityMatrix<Alloc>::iterator it_;
        
        std::vector<std::unique_ptr<imp::Matrix<Alloc>>> accDensityMatrixR_;
        std::vector<std::unique_ptr<imp::Matrix<Alloc>>> accDensityMatrixI_;
        
        
        void init(imp::EigenValues<Alloc> const& eig, double) {
            accDensityMatrixR_.resize(eig.sectorNumber() + 1);
            for(int sec = eig.sectorNumber(); sec; --sec)
                accDensityMatrixR_[sec].reset(new imp::Matrix<Alloc>(typename imp::Matrix<Alloc>::Zero(eig.at(sec).dim())));
        };
        
        void accDensityMatrix(ut::complex const sign, int sec, ut::Zahl fact, imp::Matrix<Alloc> const& matrix, imp::itf::Batcher& batcher, double) {
            imp::axpy(*accDensityMatrixR_[sec], sign.real()*fact, matrix, batcher);
        };
        
        void store(imp::EigenValues<Alloc> const& eig, jsx::value& measurements, int sec, double) {
            int const dim = eig.at(sec + 1).dim();
            std::vector<double> temp(dim*dim, .0);
            imp::axpy(temp.data(), 1., *accDensityMatrixR_[sec + 1]);
            
            int const dim0 = eig.at(sec + 1).dim0();
            std::vector<double> temp0(dim0*dim0, .0);
            for(int i = 0; i < dim; ++i)
                for(int j = 0; j < dim; ++j)
                    temp0[i*dim0 + j] = (temp[i*dim + j] + temp[j*dim + i])/2.;
            
            measurements[sec] << meas::fix(temp0, samples_);
        };
        
        void init(imp::EigenValues<Alloc> const& eig, ut::complex) {
            accDensityMatrixR_.resize(eig.sectorNumber() + 1);
            for(int sec = eig.sectorNumber(); sec; --sec)
                accDensityMatrixR_[sec].reset(new imp::Matrix<Alloc>(typename imp::Matrix<Alloc>::Zero(eig.at(sec).dim())));
            
            accDensityMatrixI_.resize(eig.sectorNumber() + 1);
            for(int sec = eig.sectorNumber(); sec; --sec)
                accDensityMatrixI_[sec].reset(new imp::Matrix<Alloc>(typename imp::Matrix<Alloc>::Zero(eig.at(sec).dim())));
        };
        
        void accDensityMatrix(ut::complex const sign, int sec, ut::Zahl fact, imp::Matrix<Alloc> const& matrix, imp::itf::Batcher& batcher, ut::complex) {
            imp::axpy(*accDensityMatrixR_[sec], sign.real()*fact, matrix, batcher);
            imp::axpy(*accDensityMatrixI_[sec], sign.imag()*fact, matrix, batcher);
        };

        void store(imp::EigenValues<Alloc> const& eig, jsx::value& measurements, int sec, ut::complex) {
            int const dim = eig.at(sec + 1).dim();
            std::vector<double> tempR(dim*dim, .0);
            imp::axpy(tempR.data(), 1., *accDensityMatrixR_[sec + 1]);
            std::vector<double> tempI(dim*dim, .0);
            imp::axpy(tempI.data(), 1., *accDensityMatrixI_[sec + 1]);
            
            int const dim0 = eig.at(sec + 1).dim0();
            std::vector<ut::complex> temp0(dim0*dim0, .0);
            for(int i = 0; i < dim; ++i)
                for(int j = 0; j < dim; ++j)
                    temp0[i*dim0 + j] = ut::complex((tempR[i*dim + j] + tempR[j*dim + i])/2., (tempI[i*dim + j] - tempI[j*dim + i])/2.);
            
            measurements[sec] << meas::fix(temp0, samples_);
        };
        
    
        void evaluate_static(ut::complex const sign, int sec, data::Data const& data, state::State& state, imp::itf::Batcher& batcher) {
            auto const& ops = imp::get<Alloc>(data.ops());
            auto const& ide = imp::get<Alloc>(data.ide());
            
            auto* bufferA = &imp::get<Alloc>(state.product()).bufferA();
            auto* bufferB = &imp::get<Alloc>(state.product()).bufferB();
            
            auto it = chain_.begin();
            
            imp::copyEvolveL(*bufferA, it->prop->at(sec), ide.mat(sec), batcher);
            while(++it + 1 != chain_.end()) {
                if(Bulla) {
                    auto const& opBulla = get<Alloc>(data.bullaOps()).at(it->flavor);
                    it->rightBulla.reset(new imp::Matrix<Alloc>(opBulla.mat(sec).I()*bufferA->J()));
                    imp::mult(*it->rightBulla, opBulla.mat(sec), *bufferA, batcher);
                }
                
                auto const& op = ops.at(it->flavor);
                if(Density) {
                    it->rightDensity.reset(new imp::Matrix<Alloc>(op.mat(sec).I()*bufferA->J()));
                    imp::mult(*it->rightDensity, op.mat(sec), *bufferA, batcher); sec = op.map(sec).sector;
                    imp::copyEvolveL(*bufferA, it->prop->at(sec), *it->rightDensity, batcher);
                } else {
                    imp::mult(*bufferB, op.mat(sec), *bufferA, batcher); sec = op.map(sec).sector;
                    imp::evolveL(it->prop->at(sec), *bufferB, batcher); std::swap(bufferA, bufferB);
                }
            }
            
            if(Bulla) {
                auto const& opBulla = get<Alloc>(data.bullaOps()).at(it->flavor);
                it->rightBulla.reset(new imp::Matrix<Alloc>(opBulla.mat(sec).I()*bufferA->J()));
                imp::mult(*it->rightBulla, opBulla.mat(sec), *bufferA, batcher);
            }
            
            ut::Zahl const fact = 1./(ut::beta()*imp::get<Alloc>(state.densityMatrix()).Z());
            
            if(Density) {
                imp::mult(*bufferB, ops.at(it->flavor).mat(sec), *bufferA, batcher); sec = ops.at(it->flavor).map(sec).sector;
                imp::density_matrix(*bufferA, ide.mat(sec), it->prop->at(sec), *bufferB, batcher);            //density matrix, last to the left:  it->prop->at(sec) *bufferB     bufferA free
                accDensityMatrix(sign, sec, fact, *bufferA, batcher, HybVal());
            } else
                sec = ops.at(it->flavor).map(sec).sector;
            
            imp::copyEvolveL(*bufferA, it->prop->at(sec), ide.mat(sec), batcher);
            do {
                if(Bulla) imp::traceAtB(nullptr, &it->acc, *bufferA, *it->rightBulla, batcher);
                
                auto const& op = ops.at(it->flavor%2 ? it->flavor - 1 : it->flavor + 1);
                imp::mult(*bufferB, op.mat(sec), *bufferA, batcher); sec = op.map(sec).sector; --it;
                
                if(Density) {
                    imp::density_matrix(*bufferA, *bufferB, it->prop->at(sec), *it->rightDensity, batcher);   //density matrix: *bufferB^T it->prop->at(sec) *it->right      bufferA free
                    accDensityMatrix(sign, sec, fact, *bufferA, batcher, HybVal());
                }
                
                evolveL(it->prop->at(sec), *bufferB, batcher); std::swap(bufferA, bufferB);
            }  while(it - 1 != chain_.begin());
            
            if(Bulla) imp::traceAtB(nullptr, &it->acc, *bufferA, *it->rightBulla, batcher);
            
            if(Density) {
                auto const& op = ops.at(it->flavor%2 ? it->flavor - 1 : it->flavor + 1);
                imp::mult(*bufferB, op.mat(sec), *bufferA, batcher); sec = op.map(sec).sector; --it;
                imp::density_matrix(*bufferA, *bufferB, it->prop->at(sec), ide.mat(sec), batcher);            //density matrix, last to the left: *bufferB^T it->prop->at(sec)
                accDensityMatrix(sign, sec, fact, *bufferA, batcher, HybVal());
            }
        };
        
        
        void evaluate_dynamic(data::Data const& data, state::State& state) {
            auto const& densityMatrix = imp::get<Alloc>(state.densityMatrix());
            auto const& dyn = *data.dyn();
            
            double N = .0;
            for(auto sec : densityMatrix)
                N += densityMatrix.weight(sec)*data.filling().at(sec);
            
            for(std::size_t i = 1; i < chain_.size(); ++i) {
                double const key = chain_[i].key; double temp = .0;
                for(std::size_t j = 1; j < chain_.size(); ++j)
                    chain_[j].flavor%2 ? temp += dyn.L(key - chain_[j].key) : temp -= dyn.L(key - chain_[j].key);
                *chain_[i].ptr += temp + (dyn.L(key) - dyn.L(key - ut::KeyMax))*N;
            }
        };
    };
}

#endif

