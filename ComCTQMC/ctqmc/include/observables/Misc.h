#ifndef OBSERVABLES_MISC_H
#define OBSERVABLES_MISC_H

#include <vector>

#include "Observable.h"
#include "../Utilities.h"
#include "../Data.h"
#include "../State.h"
#include "../../../include/measurements/Measurements.h"

namespace obs {

    struct Misc : Observable {
        Misc() = delete;
        Misc(jsx::value const& jParams, data::Data const& data) :
        print_(jParams.is("print") ? jParams("print").boolean() : false),
        accSign_(.0),
        acck_(.0),
        acckFlavor_(data.ops().flavors(), .0),
        acckHist_((jParams.is("expansion histogram") ? jParams("expansion histogram").boolean() : true) ? 1 : 0),
        accDynE_(.0),
        accDynDyn_(.0) {
        };
        Misc(Misc const&) = delete;
        Misc(Misc&&) = delete;
        Misc& operator=(Misc const&) = delete;
        Misc& operator=(Misc&&) = delete;
        ~Misc() = default;
        
        bool sample(ut::complex const sign, data::Data const& data, state::State& state, imp::itf::Batcher& batcher) {
            accSign_ += sign.real();
            
            acck_ += sign.real()*state.product().size()/2.;
            
            for(int flavor = 0; flavor < data.ops().flavors(); ++flavor)
                acckFlavor_[flavor] += sign.real()*static_cast<int>(state.config()[flavor].size());
            
            if(acckHist_.size()) {
                if(!(state.product().size()/2 < static_cast<int>(acckHist_.size())))
                    acckHist_.resize(state.product().size()/2 + 1, .0);
                acckHist_[state.product().size()/2] += sign.real();
            }
            
            accDynE_ += sign.real()*state.dyn().E();
            accDynDyn_ += sign.real()*state.dyn().dyn()*state.dyn().dyn();  //!!!!!!!!!!!!!!!!!!
            
            return true;
        };
        
        void store(data::Data const& data, jsx::value& measurements, std::int64_t samples) {
            if(print_ && samples) mpi::cout << " k: " << acck_/samples << std::endl;
            
            measurements["sign"] << meas::fix(accSign_, samples); accSign_ = .0;
            
            measurements["scalar"]["k"] << meas::fix(acck_, samples); acck_ = .0;
            
            measurements["flavor k"] << meas::fix(acckFlavor_, samples);
            for(auto& x : acckFlavor_) x = .0;
            
            if(acckHist_.size()) {
                measurements["expansion histogram"] << meas::var(acckHist_, samples);
                for(auto& x : acckHist_) x = .0;
            }
            
            measurements["dynE"] << meas::fix(accDynE_, samples); accDynE_ = .0;
            measurements["dyndyn"] << meas::fix(accDynDyn_, samples); accDynDyn_ = .0;
        };
        
        void finalize(data::Data const& data, jsx::value& measurements, std::int64_t samples) {
            store(data, measurements, samples);
        };

    private:
        bool const print_;
        
        double accSign_;
        double acck_;
        std::vector<double> acckFlavor_;
        std::vector<double> acckHist_;
        double accDynE_;
        double accDynDyn_;
    };
}

#endif
