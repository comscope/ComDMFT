#ifndef CTQMC_INCLUDE_OBSERVABLES_OBSERVABLE_H
#define CTQMC_INCLUDE_OBSERVABLES_OBSERVABLE_H

#include <vector>

#include "../Data.h"
#include "../State.h"

namespace obs {
    
    namespace itf{
        
        template<typename Value>
        struct Observable {
            virtual bool sample(Value const sign, data::Data<Value> const& data, state::State<Value>& state, jsx::value& measurements, imp::itf::Batcher<Value>& batcher) = 0;
            virtual void finalize(data::Data<Value> const& data, jsx::value& measurements) = 0;
            virtual ~Observable() = default;
        };
        
    }

}

#endif
