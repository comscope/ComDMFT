#ifndef OBSERVABLES_OBSERVABLE_H
#define OBSERVABLES_OBSERVABLE_H

#include "../Utilities.h"
#include "../Data.h"
#include "../State.h"

namespace obs {
    
    struct Observable {
        virtual bool sample(ut::complex const sign, data::Data const& data, state::State& state, imp::itf::Batcher& batcher) = 0;
        virtual void store(data::Data const& data, jsx::value& measurements, std::int64_t samples) = 0;
        virtual void finalize(data::Data const& data, jsx::value& measurements, std::int64_t samples) = 0;
        virtual ~Observable() = default;
    };
    
}

#endif
