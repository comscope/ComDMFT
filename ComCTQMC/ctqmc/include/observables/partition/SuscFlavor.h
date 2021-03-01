#ifndef CTQMC_INCLUDE_OBSERVABLES_PARTITION_SUSCFLAVOR_H
#define CTQMC_INCLUDE_OBSERVABLES_PARTITION_SUSCFLAVOR_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

#include "Exponentials.h"
#include "../Observable.h"
#include "../../Utilities.h"
#include "../../Data.h"
#include "../../State.h"
#include "../../../../include/measurements/Measurements.h"

namespace obs {
    
    namespace partition {
        
        template<typename Value>
        struct SuscFlavor : obs::itf::Observable<Value> {
            SuscFlavor() = delete;
            SuscFlavor(std::int64_t store, jsx::value const& jParams, data::Data<Value> const& data) :
            flavors_(data.ops().flavors()/2),
            nMat_(std::max(static_cast<int>(ut::beta()*jParams(cfg::partition::Worm::name())("susceptibility cutoff").real64()/(2*M_PI)), 1)),
            store_(store), samples_(0),
            acc_(flavors_, std::vector<std::vector<double>>(flavors_, std::vector<double>(nMat_, .0))),
            exponentials_(jParams, data) {
            };
            SuscFlavor(SuscFlavor const&) = delete;
            SuscFlavor(SuscFlavor&&) = delete;
            SuscFlavor& operator=(SuscFlavor const&) = delete;
            SuscFlavor& operator=(SuscFlavor&&) = delete;
            ~SuscFlavor() = default;
            
            bool sample(Value const sign, data::Data<Value> const& data, state::State<Value>& state, jsx::value& measurements, imp::itf::Batcher<Value>& batcher) {
                exponentials_.set(state.expansion());
                
                for(int f1 = 0; f1 < flavors_; ++f1)
                    for(int f2 = 0; f2 < flavors_; ++f2) {
                        auto& entry = acc_[f1][f2];
                        
                        for(auto const& key1 : state.expansion()[2*f1])
                            for(auto const& key2 : state.expansion()[2*f2]) {
                                double const time = ut::beta()*std::abs((key1 - key2)/static_cast<double>(ut::KeyMax));
                                entry[0] += ut::real(sign)*time*(ut::beta() - time);
                            }
                        
                        for(auto const& key1 : state.expansion()[2*f1 + 1])
                            for(auto const& key2 : state.expansion()[2*f2]) {
                                double const time = ut::beta()*std::abs((key1 - key2)/static_cast<double>(ut::KeyMax));
                                entry[0] -= ut::real(sign)*time*(ut::beta() - time);
                            }
                        
                        for(auto const& key1 : state.expansion()[2*f1])
                            for(auto const& key2 : state.expansion()[2*f2 + 1]) {
                                double const time = ut::beta()*std::abs((key1 - key2)/static_cast<double>(ut::KeyMax));
                                entry[0] -= ut::real(sign)*time*(ut::beta() - time);
                            }
                        
                        for(auto const& key1 : state.expansion()[2*f1 + 1])
                            for(auto const& key2 : state.expansion()[2*f2 + 1]) {
                                double const time = ut::beta()*std::abs((key1 - key2)/static_cast<double>(ut::KeyMax));
                                entry[0] += ut::real(sign)*time*(ut::beta() - time);
                            }
                        
                        auto const& exponentials1 = exponentials_.at(f1);
                        auto const& exponentials2 = exponentials_.at(f2);
                        
                        for(std::size_t n = 1; n < nMat_; ++n)
                            entry[n] += ut::real(sign*exponentials1[n]*ut::conj(exponentials2[n]));
                    }
                
                ++samples_; if(samples_%store_ == 0) store(data, measurements);
                
                return true;
            };
            
            void finalize(data::Data<Value> const& data, jsx::value& measurements) {
                store(data, measurements);
            };
            
        private:
            int const flavors_;
            std::size_t const nMat_;
            
            std::int64_t const store_;
            std::int64_t samples_;

            std::vector<std::vector<std::vector<double>>> acc_;
            Exponentials<Value> exponentials_;
            
            
            void store(data::Data<Value> const& data, jsx::value& measurements) {
                for(int f1 = 0; f1 < flavors_; ++f1)
                    for(int f2 = 0; f2 < flavors_; ++f2) {
                        measurements["susceptibility flavor"][std::to_string(f1) + "_" + std::to_string(f2)] << meas::fix(acc_[f1][f2], samples_);
                        for(auto& x : acc_[f1][f2]) x = .0;
                    }
                
                samples_ = 0;
            };
        };
        
    }
}

#endif
