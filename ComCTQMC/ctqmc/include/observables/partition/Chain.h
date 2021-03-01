#ifndef CTQMC_INCLUDE_OBSERVABLES_PARTITION_CHAIN_H
#define CTQMC_INCLUDE_OBSERVABLES_PARTITION_CHAIN_H

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include "DensityMatrix.h"
#include "../Observable.h"
#include "../../Utilities.h"
#include "../../Data.h"
#include "../../State.h"
#include "../../../../include/measurements/Measurements.h"

namespace obs {
    
    namespace partition {
        
        // Todo: search for conjugate as a function of the pointer because flavor will disapear form node
        
        template<typename Mode, typename Value, bool Density, bool Bulla>
        struct Chain : obs::itf::Observable<Value> {
            Chain() = delete;
            Chain(std::int64_t, jsx::value const& jParams, data::Data<Value> const& data) :
            phase_(Phase::Prepare), samples_(0) {
                if(Density) densityMatrix_.reset(new DensityMatrix<Mode, Value>(jParams, data));
            };
            Chain(Chain const&) = delete;
            Chain(Chain&&) = delete;
            Chain operator=(Chain const&) = delete;
            Chain& operator=(Chain&&) = delete;
            ~Chain() = default;
            
            bool sample(Value const sign, data::Data<Value> const& data, state::State<Value>& state, jsx::value& measurements, imp::itf::Batcher<Value>& batcher) { // scheisse das ...
                auto const& densityMatrix = imp::get<Mode>(state.densityMatrix());
                auto& product = imp::get<Mode>(state.product());
                
                if(product.size()) {
                    switch(phase_) {
                        case Phase::Prepare:
                            prepare_chain(data, state);
                            
                            it_ = densityMatrix.begin();
                            
                            phase_ = Phase::Calculate;
                            
                        case Phase::Calculate:
                            if(it_ != densityMatrix.end()) {
                                evaluate_sec(sign, *it_, data, state, batcher);
                                ++it_; return false;
                            }
                            
                            if(Bulla) {
                                for(auto it = chain_.begin() + 1; it != chain_.end(); ++it) *it->ptr = (it->acc/densityMatrix.Z()).get();
                                if(data.dyn() != nullptr) dynamic(data, state);
                            }
                            
                            chain_.clear();
                            
                            phase_ = Phase::Prepare;
                    }
                } else if(Density) {
                    auto const& ide = imp::get<Mode>(data.ide());
                    auto const& eig = imp::get<Mode>(data.eig());
                    
                    for(auto const sec : densityMatrix) {
                        imp::density_matrix(product.bufferA(), ide.mat(sec), product.first()->prop()->at(sec), ide.mat(sec), eig.at(sec), batcher);
                        imp::add(densityMatrix_->mat(sec), (sign/ut::beta())/densityMatrix.Z(), product.bufferA(), batcher);
                    }
                }
                
                ++samples_; return true;
            };
            
            
            void finalize(data::Data<Value> const& data, jsx::value& measurements) {
                if(Density) densityMatrix_->store(data, measurements["density matrix"], samples_);
            };
            
        private:
            
            enum class Phase { Prepare, Calculate };
            
            struct Entry {
                Entry(ut::KeyType key, int flavor, imp::Propagator<Mode> const* prop, Value* ptr) :
                key(key), flavor(flavor), prop(prop), acc(.0), ptr(ptr) {
                };
                
                ut::KeyType const key;
                int const flavor;
                imp::Propagator<Mode> const* const prop;
                
                ut::Zahl<Value> acc;
                Value* const ptr;
                
                std::unique_ptr<imp::Matrix<Mode, Value>> rightDensity, rightBulla;
            };
            
            
            Phase phase_;
            std::int64_t samples_;
            
            std::vector<Entry> chain_;
            typename imp::DensityMatrix<Mode, Value>::iterator it_;
            
            std::unique_ptr<DensityMatrix<Mode, Value>> densityMatrix_;
            
            
            void prepare_chain(data::Data<Value> const& data, state::State<Value>& state) {
                auto const& product = imp::get<Mode>(state.product());
                
                if(Bulla) {
                    std::map<ut::KeyType, Value*> ptrs;  ptrs[0] = nullptr;
                    for(auto const& bath : state.baths()) {
                        for(auto const& opL : bath.opsL()) ptrs[opL.key()] = &opL.bulla();
                        for(auto const& opR : bath.opsR()) ptrs[opR.key()] = &opR.bulla();
                    }
                    
                    auto itPtr = ptrs.begin();
                    for(auto it = product.first(); it != product.last(); it = it.next(0), ++itPtr)
                        chain_.emplace_back(it.key(), it->flavor, it->prop(), itPtr->second);
                } else {
                    for(auto it = product.first(); it != product.last(); it = it.next(0))
                        chain_.emplace_back(it.key(), it->flavor, it->prop(), nullptr);
                }
            }
            
            
            void evaluate_sec(Value const sign, int sec, data::Data<Value> const& data, state::State<Value>& state, imp::itf::Batcher<Value>& batcher) {
                auto const& ops = imp::get<Mode>(data.ops());
                auto const& ide = imp::get<Mode>(data.ide());
                auto const& eig = imp::get<Mode>(data.eig());
                auto const& bullaOps = imp::get<Mode>(data.template opt<imp::itf::BullaOperators<Value>>());
                
                auto* bufferA = &imp::get<Mode>(state.product()).bufferA();
                auto* bufferB = &imp::get<Mode>(state.product()).bufferB();
                
                auto it = chain_.begin();
                
                imp::copyEvolveL(*bufferA, it->prop->at(sec), ide.mat(sec), batcher);
                while(++it + 1 != chain_.end()) {
                    if(Bulla) {
                        if(bullaOps.at(it->flavor).map(sec).sector != 0) {
                            it->rightBulla.reset(new imp::Matrix<Mode, Value>(bullaOps.at(it->flavor).mat(sec).I()*bufferA->J()));
                            imp::mult(*it->rightBulla, bullaOps.at(it->flavor).mat(sec), *bufferA, batcher);
                        } else
                            it->rightBulla.reset(nullptr);
                    }
                    
                    auto const& op = imp::get<Mode, Value>(ops.at(it->flavor));
                    if(Density) {
                        it->rightDensity.reset(new imp::Matrix<Mode, Value>(op.mat(sec).I()*bufferA->J()));
                        imp::mult(*it->rightDensity, op.mat(sec), *bufferA, batcher); sec = op.map(sec).sector;
                        imp::copyEvolveL(*bufferA, it->prop->at(sec), *it->rightDensity, batcher);
                    } else {
                        imp::mult(*bufferB, op.mat(sec), *bufferA, batcher); sec = op.map(sec).sector;
                        imp::evolveL(it->prop->at(sec), *bufferB, batcher); std::swap(bufferA, bufferB);
                    }
                }
                
                if(Bulla) {
                    if(bullaOps.at(it->flavor).map(sec).sector != 0) {
                        it->rightBulla.reset(new imp::Matrix<Mode, Value>(bullaOps.at(it->flavor).mat(sec).I()*bufferA->J()));
                        imp::mult(*it->rightBulla, bullaOps.at(it->flavor).mat(sec), *bufferA, batcher);
                    } else
                        it->rightBulla.reset(nullptr);
                }
                
                auto const fact = (sign/ut::beta())/imp::get<Mode>(state.densityMatrix()).Z();
                
                if(Density) {
                    auto const& op = imp::get<Mode, Value>(ops.at(it->flavor));
                    imp::mult(*bufferB, op.mat(sec), *bufferA, batcher); sec = op.map(sec).sector;
                    imp::density_matrix(*bufferA, ide.mat(sec), it->prop->at(sec), *bufferB, eig.at(sec), batcher);            //density matrix, last to the left:  it->prop->at(sec) *bufferB     bufferA free
                    imp::add(densityMatrix_->mat(sec), fact, *bufferA, batcher);
                } else {
                    auto const& op = imp::get<Mode, Value>(ops.at(it->flavor));
                    sec = op.map(sec).sector;
                }
                
                imp::copyEvolveL(*bufferA, it->prop->at(sec), ide.mat(sec), batcher);
                do {
                    if(Bulla && it->rightBulla.get() != nullptr) imp::traceAtB(static_cast<ut::Zahl<Value>*>(nullptr), &it->acc, *bufferA, *it->rightBulla, batcher);
                    
                    auto const& op = imp::get<Mode, Value>(ops.at(it->flavor%2 ? it->flavor - 1 : it->flavor + 1));
                    imp::mult(*bufferB, op.mat(sec), *bufferA, batcher); sec = op.map(sec).sector; --it;
                    
                    if(Density) {
                        imp::density_matrix(*bufferA, *bufferB, it->prop->at(sec), *it->rightDensity, eig.at(sec), batcher);   //density matrix: *bufferB^t it->prop->at(sec) *it->right      bufferA free
                        imp::add(densityMatrix_->mat(sec), fact, *bufferA, batcher);
                    }
                    
                    evolveL(it->prop->at(sec), *bufferB, batcher); std::swap(bufferA, bufferB);
                }  while(it - 1 != chain_.begin());
                
                if(Bulla && it->rightBulla.get() != nullptr) imp::traceAtB(static_cast<ut::Zahl<Value>*>(nullptr), &it->acc, *bufferA, *it->rightBulla, batcher);
                
                if(Density) {
                    auto const& op = imp::get<Mode, Value>(ops.at(it->flavor%2 ? it->flavor - 1 : it->flavor + 1));
                    imp::mult(*bufferB, op.mat(sec), *bufferA, batcher); sec = op.map(sec).sector; --it;
                    imp::density_matrix(*bufferA, *bufferB, it->prop->at(sec), ide.mat(sec), eig.at(sec), batcher);            //density matrix, last to the left: *bufferB^t it->prop->at(sec)
                    imp::add(densityMatrix_->mat(sec), fact, *bufferA, batcher);
                }
            };
            
            void dynamic(data::Data<Value> const& data, state::State<Value>& state) {
                auto const& dynFunc = *data.dyn();
                auto const& densityMatrix = imp::get<Mode>(state.densityMatrix());
                auto const& dynWeight = state.dyn();
                
                std::vector<Value> boundary(data.ops().flavors(), .0);
                
                for(auto sec : densityMatrix)
                    for(int flavor = 0; flavor < data.ops().flavors(); ++flavor)
                        boundary[flavor] += densityMatrix.weight(sec)*dynFunc.D0Qf(sec, flavor);
                
                for(std::size_t i = 1; i < chain_.size(); ++i)
                    *chain_[i].ptr += boundary[chain_[i].flavor] - dynWeight.fkinks(chain_[i].flavor, chain_[i].key);
            };
        };
        
    }
}

#endif

