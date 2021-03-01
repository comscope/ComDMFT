#ifndef CTQMC_INCLUDE_OBSERVABLES_WORM_SETUP_H
#define CTQMC_INCLUDE_OBSERVABLES_WORM_SETUP_H

#include "Observable.h"

#include "../Observables.h"
#include "../../config/Worms.h"

namespace obs {

    template<typename Mode, typename Worm, worm::MeasType measType, typename Value>
    void add_obs(jsx::value const& jWorm, data::Data<Value>& data, WormObservables<Value>& observables)
    {
        std::int64_t const sweep = jWorm.is("sweep") ? jWorm("sweep").int64() : 50;
        std::int64_t const store = jWorm.is("store") ? jWorm("store").int64() : 100;
        
        if(jWorm.is("basis") ? jWorm("basis").string() == "matsubara" : true)
            observables.template add<worm::Observable<Mode, Value, worm::FuncType::Matsubara, measType, Worm>>(sweep, store, jWorm, data);
        else if(jWorm.is("basis") ? jWorm("basis").string() == "legendre" : true)
            observables.template add<worm::Observable<Mode, Value, worm::FuncType::Legendre, measType, Worm>>(sweep, store, jWorm, data);
        else
            throw std::runtime_error("Unknown basis option");
    }

    
    template<typename Mode, typename Value, typename Worm>
    void setup_worm_obs(jsx::value const& jParams, data::Data<Value>& data, std::unique_ptr<WormObservables<Value>>& observables, ut::wrap<Worm>)
    {
        auto const& jWorm = jParams(Worm::name());
        
        observables.reset(new WormObservables<Value>(Worm::name()));
        
        if(jWorm.is("static"))
            add_obs<Mode, Worm, worm::MeasType::Static>(jWorm, data, *observables);
        
        if(jWorm.is("dynamic"))
            add_obs<Mode, Worm, worm::MeasType::Dynamic>(jWorm, data, *observables);
    }
    
}

#endif
