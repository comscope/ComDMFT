#ifndef CTQMC_INCLUDE_UPDATES_EXPANSION_SETUP_H
#define CTQMC_INCLUDE_UPDATES_EXPANSION_SETUP_H

#include <vector>
#include <array>

#include "InsertErasePairCsc.h"
#include "InsertErasePairStd.h"
#include "InsertEraseQuadStd.h"

#include "../../markovchain/Update.h"
#include "../../markovchain/MarkovChain.h"

// Scheisse kack 4 operator updates immer no nit implementiert Himmelhergottsakramentzifixhallelujaleckstmiamarscheisseglumpvereckts !

namespace upd {
    
    namespace expansion {
        
        template<typename Space, typename Mode, typename Value>
        void setup_updates(jsx::value const& jParams, data::Data<Value> const& data, state::State<Value> const& state, mch::MarkovChain<Value>& markovChain) {
            
            mpi::cout << "Setting expansion updates for " + Space::name() + " space ... ";
            
            markovChain.add(mch::unique_update_ptr<Value>(new upd::Generic< InsertCsc<Space>, Mode, Value >(1., jParams, data)),
                            mch::unique_update_ptr<Value>(new upd::Generic< EraseCsc<Space>, Mode, Value >(1., jParams, data)));
            
            if (jParams.is("quad insert") and jParams("quad insert").boolean()){
                        markovChain.add(mch::unique_update_ptr<Value>(new upd::Generic< QuadInsertStd<Space>, Mode, Value >(1., jParams, data)),
                                        mch::unique_update_ptr<Value>(new upd::Generic< QuadEraseStd<Space>, Mode, Value >(1., jParams, data)));
            }
            
            mpi::cout << "Ok" << std::endl;
            
        };
        
    }
    
}

#endif
