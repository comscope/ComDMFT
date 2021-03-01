#ifndef CTQMC_INCLUDE_OBSERVABLES_PARTITION_MISC_H
#define CTQMC_INCLUDE_OBSERVABLES_PARTITION_MISC_H

#include <vector>

#include "../Observable.h"
#include "../../Utilities.h"
#include "../../Data.h"
#include "../../State.h"
#include "../../../../include/measurements/Measurements.h"

namespace obs {
    
    namespace partition {
        
        template<typename Value>
        struct Misc : obs::itf::Observable<Value> {
            Misc() = delete;
            Misc(std::int64_t store, jsx::value const& jParams, data::Data<Value> const& data) :
            print_(jParams.is("print") ? jParams("print").boolean() : false),
            store_(store), samples_(0),
            accSign_(.0),
            acck_(.0),
            acckFlavor_(data.ops().flavors(), .0),
            acckHist_((jParams(cfg::partition::Worm::name()).is("expansion histogram") ? jParams(cfg::partition::Worm::name())("expansion histogram").boolean() : true) ? 1 : 0),
            accDynEnergy_(.0) {
            };
            Misc(Misc const&) = delete;
            Misc(Misc&&) = delete;
            Misc& operator=(Misc const&) = delete;
            Misc& operator=(Misc&&) = delete;
            ~Misc() = default;
            
            bool sample(Value const sign, data::Data<Value> const& data, state::State<Value>& state, jsx::value& measurements, imp::itf::Batcher<Value>& batcher) {
                accSign_ += ut::real(sign);
                
                acck_ += ut::real(sign)*state.product().size()/2.;
                
                for(int flavor = 0; flavor < data.ops().flavors(); ++flavor)
                    acckFlavor_[flavor] += ut::real(sign)*static_cast<int>(state.expansion()[flavor].size());
                
                if(acckHist_.size()) {
                    if(!(state.product().size()/2 < static_cast<int>(acckHist_.size())))
                        acckHist_.resize(state.product().size()/2 + 1, .0);
                    acckHist_[state.product().size()/2] += ut::real(sign);
                }
                
                accDynEnergy_ += ut::real(sign)*state.dyn().energy();
                    
                ++samples_; if(samples_%store_ == 0) store(data, measurements);
                
                return true;
            };
            
            void finalize(data::Data<Value> const& data, jsx::value& measurements) {
                store(data, measurements);
            };
            
        private:
            bool const print_;
            
            std::int64_t const store_;
            std::int64_t samples_;
            
            double accSign_;
            double acck_;
            std::vector<double> acckFlavor_;
            std::vector<double> acckHist_;
            double accDynEnergy_;
            
            
            void store(data::Data<Value> const& data, jsx::value& measurements) {
                if(print_ && samples_) mpi::cout << " k: " << acck_/samples_ << std::endl;
                
                measurements["sign"] << meas::fix(accSign_, samples_); accSign_ = .0;
                
                measurements["scalar"]["k"] << meas::fix(acck_, samples_); acck_ = .0;
                
                measurements["flavor k"] << meas::fix(acckFlavor_, samples_);
                for(auto& x : acckFlavor_) x = .0;
                
                if(acckHist_.size()) {
                    measurements["expansion histogram"] << meas::var(acckHist_, samples_);
                    for(auto& x : acckHist_) x = .0;
                }
                
                measurements["dyn energy"] << meas::fix(accDynEnergy_, samples_); accDynEnergy_ = .0;
                
                samples_ = 0;
            };
        };
        
    }
}

#endif
