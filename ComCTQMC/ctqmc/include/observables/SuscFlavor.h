#ifndef OBSERVABLES_SUSCFLAVOR_H
#define OBSERVABLES_SUSCFLAVOR_H

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
    
    struct SuscFlavor : Observable {
        SuscFlavor() = delete;
        SuscFlavor(jsx::value const& jParams, data::Data const& data) :
        flavors_(data.ops().flavors()/2),
        nMat_(std::max(static_cast<int>(ut::beta()*jParams("susceptibility cutoff").real64()/(2*M_PI)), 1)),
        acc_(flavors_, std::vector<std::vector<double>>(flavors_, std::vector<double>(nMat_, .0))) {
        };
        SuscFlavor(SuscFlavor const&) = delete;
        SuscFlavor(SuscFlavor&&) = delete;
        SuscFlavor& operator=(SuscFlavor const&) = delete;
        SuscFlavor& operator=(SuscFlavor&&) = delete;
        ~SuscFlavor() = default;
        
        bool sample(ut::complex const sign, data::Data const& data, state::State& state, imp::itf::Batcher& batcher) {
            for(int f1 = 0; f1 < flavors_; ++f1)
                for(int f2 = 0; f2 < flavors_; ++f2) {
                    auto& entry = acc_[f1][f2];
                    
                    for(auto const& c1 : state.config()[2*f1])
                        for(auto const& c2 : state.config()[2*f2]) {
                            double const time = ut::beta()*std::abs((c1.key() - c2.key())/static_cast<double>(ut::KeyMax));
                            entry[0] += sign.real()*time*(ut::beta() - time);
                        }
                    
                    for(auto const& c1 : state.config()[2*f1 + 1])
                        for(auto const& c2 : state.config()[2*f2]) {
                            double const time = ut::beta()*std::abs((c1.key() - c2.key())/static_cast<double>(ut::KeyMax));
                            entry[0] -= sign.real()*time*(ut::beta() - time);
                        }
                    
                    for(auto const& c1 : state.config()[2*f1])
                        for(auto const& c2 : state.config()[2*f2 + 1]) {
                            double const time = ut::beta()*std::abs((c1.key() - c2.key())/static_cast<double>(ut::KeyMax));
                            entry[0] -= sign.real()*time*(ut::beta() - time);
                        }
                    
                    for(auto const& c1 : state.config()[2*f1 + 1])
                        for(auto const& c2 : state.config()[2*f2 + 1]) {
                            double const time = ut::beta()*std::abs((c1.key() - c2.key())/static_cast<double>(ut::KeyMax));
                            entry[0] += sign.real()*time*(ut::beta() - time);
                        }

                    auto const& exponentials1 = state.exponentials().at(f1);
                    auto const& exponentials2 = state.exponentials().at(f2);
                    
                    for(std::size_t n = 1; n < nMat_; ++n)
                        entry[n] += (sign*exponentials1[n]*std::conj(exponentials2[n])).real();
                }

            return true;
        };
        
        void store(data::Data const& data, jsx::value& measurements, std::int64_t samples) {
            for(int f1 = 0; f1 < flavors_; ++f1)
                for(int f2 = 0; f2 < flavors_; ++f2) {
                    measurements["susceptibility flavor"][std::to_string(f1) + "_" + std::to_string(f2)] << meas::fix(acc_[f1][f2], samples);
                    for(auto& x : acc_[f1][f2]) x = .0;
                }
        };
        
        void finalize(data::Data const& data, jsx::value& measurements, std::int64_t samples) {
            store(data, measurements, samples);
        };

    private:
        int const flavors_;
        std::size_t const nMat_;

        std::vector<std::vector<std::vector<double>>> acc_;
    };

}

#endif 
