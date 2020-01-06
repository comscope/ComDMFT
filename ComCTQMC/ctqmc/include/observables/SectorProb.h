#ifndef OBSERVABLES_SECTORPROB_H
#define OBSERVABLES_SECTORPROB_H

#include <vector>

#include "Observable.h"
#include "../Utilities.h"
#include "../Data.h"
#include "../State.h"
#include "../../../include/measurements/Measurements.h"

namespace obs {
    
    template<typename Alloc>
    struct SectorProb  : Observable {
        SectorProb() = delete;
        SectorProb(jsx::value const& jParams, data::Data const& data) :
        acc_(data.eig().sectorNumber() + 1, .0) {
        };
        SectorProb(SectorProb const&) = delete;
        SectorProb(SectorProb&&) = delete;
        SectorProb& operator=(SectorProb const&) = delete;
        SectorProb& operator=(SectorProb&&) = delete;
        ~SectorProb() = default;
        
        bool sample(ut::complex const sign, data::Data const& data, state::State& state, imp::itf::Batcher& batcher) {
            auto& product = imp::get<Alloc>(state.product());
            auto& densityMatrix = imp::get<Alloc>(state.densityMatrix());
            
            std::vector<std::pair<int, int>> sector;  for(auto s : densityMatrix) sector.push_back(std::make_pair(s, s));

            double prev = .0;
            for(auto n = product.first().next(0); n != product.last(); n = n.next(0)) {
                double const present = n.key()/static_cast<double>(ut::KeyMax);
                
                double const diff = present - prev;
                for(auto& s : sector) {
                    acc_[s.second] += sign.real()*densityMatrix.weight(s.first)*diff;
                    s.second = n->op0->map(s.second).sector;
                }
                prev = present;
            }

            double diff = 1. - prev;
            for(auto s : sector)
                acc_[s.second] += sign.real()*densityMatrix.weight(s.first)*diff;

            return true;
        };
        
        void store(data::Data const& data, jsx::value& measurements, std::int64_t samples) {
            acc_.erase(acc_.begin());
            measurements["sector prob"] << meas::fix(acc_, samples);  acc_.clear(); 
            acc_.resize(data.eig().sectorNumber() + 1, .0);
        };
        
        void finalize(data::Data const& data, jsx::value& measurements, std::int64_t samples) {
            store(data, measurements, samples);
        };
        
    private:
        std::vector<double> acc_;
    };
}

#endif
