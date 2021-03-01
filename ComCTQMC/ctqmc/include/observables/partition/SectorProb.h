#ifndef CTQMC_INCLUDE_OBSERVABLES_PARTITION_SECTORPROB_H
#define CTQMC_INCLUDE_OBSERVABLES_PARTITION_SECTORPROB_H

#include <vector>

#include "../Observable.h"
#include "../../Utilities.h"
#include "../../Data.h"
#include "../../State.h"
#include "../../../../include/measurements/Measurements.h"

namespace obs {
    
    namespace partition {
        
        template<typename Mode, typename Value>
        struct SectorProb  : obs::itf::Observable<Value> {
            SectorProb() = delete;
            SectorProb(std::int64_t store, jsx::value const& jParams, data::Data<Value> const& data) :
            store_(store), samples_(0),
            acc_(data.eig().sectorNumber() + 1, .0) {
            };
            SectorProb(SectorProb const&) = delete;
            SectorProb(SectorProb&&) = delete;
            SectorProb& operator=(SectorProb const&) = delete;
            SectorProb& operator=(SectorProb&&) = delete;
            ~SectorProb() = default;
            
            bool sample(Value const sign, data::Data<Value> const& data, state::State<Value>& state, jsx::value& measurements, imp::itf::Batcher<Value>& batcher) {
                auto& product = imp::get<Mode>(state.product());
                
                std::vector<std::pair<int, int>> sector;  for(auto s : state.densityMatrix()) sector.push_back(std::make_pair(s, s));
                
                double prev = .0;
                for(auto n = product.first().next(0); n != product.last(); n = n.next(0)) {
                    double const present = n.key()/static_cast<double>(ut::KeyMax);
                    
                    double const diff = present - prev;
                    for(auto& s : sector) {
                        acc_[s.second] += ut::real(sign*state.densityMatrix().weight(s.first))*diff;
                        s.second = n->op0->map(s.second).sector;
                    }
                    prev = present;
                }
                
                double diff = 1. - prev;
                for(auto s : sector)
                    acc_[s.second] += ut::real(sign*state.densityMatrix().weight(s.first))*diff;
                
                ++samples_; if(samples_%store_ == 0) store(data, measurements);
                
                return true;
            };
            
            void finalize(data::Data<Value> const& data, jsx::value& measurements) {
                store(data, measurements);
            };
            
        private:
            std::int64_t const store_;
            std::int64_t samples_;
            
            std::vector<double> acc_;
            
            
            void store(data::Data<Value> const& data, jsx::value& measurements) {
                acc_.erase(acc_.begin());
                measurements["sector prob"] << meas::fix(acc_, samples_);  acc_.clear();
                acc_.resize(data.eig().sectorNumber() + 1, .0);
                
                samples_ = 0;
            };
        };
        
    }
}

#endif
