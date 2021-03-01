#ifndef CTQMC_INCLUDE_OBSERVABLES_PARTITION_DENSITYMATRIX_H
#define CTQMC_INCLUDE_OBSERVABLES_PARTITION_DENSITYMATRIX_H

#include <vector>

#include "../Observable.h"
#include "../../Utilities.h"
#include "../../Data.h"
#include "../../State.h"
#include "../../../../include/measurements/Measurements.h"

namespace obs {
    
    namespace partition {
        
        template<typename Mode, typename Value>
        struct DensityMatrix {
            DensityMatrix() = delete;
            DensityMatrix(jsx::value const& jParams, data::Data<Value> const& data) {
                auto const& eig = imp::get<Mode>(data.eig());
                
                data_.resize(eig.sectorNumber() + 1);
                for(int sec = eig.sectorNumber(); sec; --sec)
                    data_[sec].reset(new imp::Matrix<Mode, Value>(typename imp::Matrix<Mode, Value>::Zero(eig.at(sec).dim())));
            };
            DensityMatrix(DensityMatrix const&) = delete;
            DensityMatrix(DensityMatrix&&) = delete;
            DensityMatrix& operator=(DensityMatrix const&) = delete;
            DensityMatrix& operator=(DensityMatrix&&) = delete;
            ~DensityMatrix() = default;
            
            imp::Matrix<Mode, Value>& mat(int sec) {
                return *data_[sec];
            };
            
            void store(data::Data<Value> const& data, jsx::value& measurements, std::int64_t samples) {
                if(!measurements.is<jsx::array_t>())
                    measurements = jsx::array_t(data.eig().sectorNumber());
                
                auto const& eig = imp::get<Mode>(data.eig());
                
                for(int sec = 0; sec < eig.sectorNumber(); ++sec) {
                    int const dim = eig.at(sec + 1).dim();
                    int const dim0 = eig.at(sec + 1).dim0();
                    
                    std::vector<Value> temp(dim*dim, .0);
                    imp::add(temp.data(), Value{1.}, *data_[sec + 1]);

                    std::vector<Value> temp0(dim0*dim0, .0);
                    for(int i = 0; i < dim; ++i)
                        for(int j = 0; j < dim; ++j)
                            temp0[i*dim0 + j] = (temp[i*dim + j] + ut::conj(temp[j*dim + i]))/2.;   // row or column major ??????
                    
                    measurements[sec] << meas::fix(temp0, samples);
                }
            };
            
        private:
            
            std::vector<std::unique_ptr<imp::Matrix<Mode, Value>>> data_;
        };
        
        
        template<typename Mode, typename Value>
        struct DensityMatrixStatic : obs::itf::Observable<Value> {
            DensityMatrixStatic() = delete;
            DensityMatrixStatic(std::int64_t, jsx::value const& jParams, data::Data<Value> const& data) :
            phase_(Phase::Sample), samples_(0),
            accDensityMatrix_(jParams, data) {
            };
            DensityMatrixStatic(DensityMatrixStatic const&) = delete;
            DensityMatrixStatic(DensityMatrixStatic&&) = delete;
            DensityMatrixStatic& operator=(DensityMatrixStatic const&) = delete;
            DensityMatrixStatic& operator=(DensityMatrixStatic&&) = delete;
            ~DensityMatrixStatic() = default;
            
            bool sample(Value const sign, data::Data<Value> const& data, state::State<Value>& state, jsx::value& measurements, imp::itf::Batcher<Value>& batcher) {
                auto const& densityMatrix = imp::get<Mode>(state.densityMatrix());
                
                switch(phase_) {
                    case Phase::Sample:
                        for(auto sec : densityMatrix)
                            imp::add(accDensityMatrix_.mat(sec), sign/densityMatrix.Z(), densityMatrix.mat(sec), batcher);
                        
                        phase_ = Phase::Finalize; return false;
                        
                    case Phase::Finalize:
                        phase_ = Phase::Sample;
                }
                
                ++samples_; return true;
            }
            
            void finalize(data::Data<Value> const& data, jsx::value& measurements) {
                accDensityMatrix_.store(data, measurements["density matrix"], samples_);
            };
            
        private:
            enum class Phase { Sample, Finalize };
            
            Phase phase_;
            std::int64_t samples_;
            DensityMatrix<Mode, Value> accDensityMatrix_;
            
        };
        
        
        template<typename Mode, typename Value>
        struct DensityMatrixDynamic : obs::itf::Observable<Value> {
            DensityMatrixDynamic() = delete;
            DensityMatrixDynamic(std::int64_t store, jsx::value const& jParams, data::Data<Value> const& data) :
            store_(store), phase_(Phase::Sample), samples_(0), samples0_(0),
            accQQ_(data.dyn()->size(), std::vector<double>(data.dyn()->size(), .0)),
            accDensityMatrices_(data.dyn()->size()) {
                for(auto& entry : accDensityMatrices_) entry.reset(new DensityMatrix<Mode, Value>(jParams, data));
            };
            DensityMatrixDynamic(DensityMatrixDynamic const&) = delete;
            DensityMatrixDynamic(DensityMatrixDynamic&&) = delete;
            DensityMatrixDynamic& operator=(DensityMatrixDynamic const&) = delete;
            DensityMatrixDynamic& operator=(DensityMatrixDynamic&&) = delete;
            ~DensityMatrixDynamic() = default;
            
            bool sample(Value const sign, data::Data<Value> const& data, state::State<Value>& state, jsx::value& measurements, imp::itf::Batcher<Value>& batcher) {
                auto const& dynFunc = *data.dyn();
                auto const& densityMatrix = imp::get<Mode>(state.densityMatrix());
                auto const& dynWeight = state.dyn();
                
                switch(phase_) {
                    case Phase::Sample:
                        for(int qn = 0; qn < data.dyn()->size(); ++qn)
                            for(auto sec : densityMatrix) {
                                auto const fact = (dynFunc.D0Qq(sec, qn) - dynWeight.qkinks(qn))*sign/densityMatrix.Z();
                                imp::add(accDensityMatrices_[qn]->mat(sec), fact, densityMatrix.mat(sec), batcher);
                            }
                        
                        for(int qnI = 0; qnI < data.dyn()->size(); ++qnI)
                            for(int qnJ = 0; qnJ < data.dyn()->size(); ++qnJ)
                                for(auto sec : densityMatrix) {
                                    auto const factI = (dynFunc.D0Qq(sec, qnI) - dynWeight.qkinks(qnI));
                                    auto const factJ = (dynFunc.D0Qq(sec, qnJ) - dynWeight.qkinks(qnJ));
                                    accQQ_[qnI][qnJ] += ut::real(sign*densityMatrix.weight(sec))*factI*factJ;
                                }

                        phase_ = Phase::Finalize; return false;
                        
                    case Phase::Finalize:
                        
                        phase_ = Phase::Sample;
                }
                
                ++samples_; if(samples_%store_ == 0) store(data, measurements);
                
                return true;
            }
            
            void finalize(data::Data<Value> const& data, jsx::value& measurements) {
                store(data, measurements);
                
                measurements["density matrix dyn"] = jsx::array_t(data.dyn()->size());
                for(int qn = 0; qn < data.dyn()->size(); ++qn)
                    accDensityMatrices_[qn]->store(data, measurements["density matrix dyn"][qn], samples0_);
            };
            
        private:
            enum class Phase { Sample, Finalize };
            
            std::int64_t const store_;
            
            Phase phase_;
            std::int64_t samples_, samples0_;
            std::vector<std::vector<double>> accQQ_;
            std::vector<std::unique_ptr<DensityMatrix<Mode, Value>>> accDensityMatrices_;
            
            
            void store(data::Data<Value> const& data, jsx::value& measurements) {
                samples0_ += samples_;
                
                if(!measurements["QQ"].is<jsx::array_t>()) measurements["QQ"] = jsx::array_t(data.dyn()->size());
                
                for(int qn = 0; qn < data.dyn()->size(); ++qn) {
                    measurements["QQ"][qn] << meas::fix(accQQ_[qn], samples_);
                    std::fill(accQQ_[qn].begin(), accQQ_[qn].end(), .0);
                }

                samples_ = .0;
            }
            
        };
  
    }
}

#endif
