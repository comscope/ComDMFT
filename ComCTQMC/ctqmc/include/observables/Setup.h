#ifndef CTQMC_INCLUDE_OBSERVABLES_SETUP_H
#define CTQMC_INCLUDE_OBSERVABLES_SETUP_H

#include <vector>

#include "partition/Setup.h"
#include "worm/Setup.h"


namespace obs {

    
    template<typename Mode, typename Value>
    struct setup_worm_obs_functor {
        template<typename W> using get_index = cfg::get_index<W, cfg::Worm>;
        
        template<typename W>
        void operator()(ut::wrap<W> w, jsx::value const& jParams, data::Data<Value>& data, Observables<Value>& observables) const {
            
            if(jParams.is(W::name())) {
                mpi::cout << "Begin setting up " + W::name() + " observables" << std::endl;
                
                setup_worm_obs<Mode>(jParams, data, observables[get_index<W>::value], w);
                
                mpi::cout << "End setting up " + W::name() + " observables" << std::endl;
            }
            
        }
    };
    

    template<typename Mode, typename Value>
    void setup_obs(jsx::value const& jParams, data::Data<Value>& data, Observables<Value>& observables) {

        cfg::for_each_type<cfg::Worm>::apply(setup_worm_obs_functor<Mode, Value>(), jParams, data, observables);
        
    }

}

#endif
