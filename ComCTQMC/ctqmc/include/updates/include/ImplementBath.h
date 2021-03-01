#ifndef CTQMC_INCLUDE_UPDATES_INCLUDE_IMPLEMENTBATH_H
#define CTQMC_INCLUDE_UPDATES_INCLUDE_IMPLEMENTBATH_H


#include "../../Data.h"
#include "../../State.h"


namespace upd {
    
    
    template<typename Update, typename Value, typename = void>
    struct ImplementBath {
        
        double ratio(Update const& update, data::Data<Value> const& data, state::State<Value>& state) {
            return 1.;
        };
        
        int accept(Update const& update, data::Data<Value> const& data, state::State<Value>& state) {
            return 1;
        };
        
        void reject(Update const& update, data::Data<Value> const& data, state::State<Value>& state) {
        };
        
    };
    
    
    template<typename Update, typename Value>
    struct ImplementBath<Update, Value, ut::void_t<decltype(&Update::template bath<Value>)>> {
        
        double ratio(Update const& update, data::Data<Value> const& data, state::State<Value>& state) {
            update.bath(data, state.baths()); double temp = 1.;
            for(auto& bath : state.baths()) temp *= bath.ratio(data.hyb());
            return temp;
        }
        
        void accept(Update const& update, data::Data<Value> const& data, state::State<Value>& state) {
            for(auto& bath : state.baths()) state.signTimeOrder() *= bath.accept(data.hyb());
        }
        
        void reject(Update const& update, data::Data<Value> const& data, state::State<Value>& state) {
            for(auto& bath : state.baths()) bath.reject(data.hyb());
        }
        
    };

}

#endif
