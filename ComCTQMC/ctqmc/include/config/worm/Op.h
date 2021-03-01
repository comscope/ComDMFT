#ifndef CTQMC_INCLUDE_CONFIG_WORM_OP_H
#define CTQMC_INCLUDE_CONFIG_WORM_OP_H


#include "Basic.h"
#include "../../Utilities.h"
#include "../../../../include/JsonX.h"


namespace cfg {

    
    template<std::size_t type>
    struct Op : Flavor, FermionicTime {
        Op() = default;
        Op(int flavor, ut::KeyType key = 0) : Flavor(flavor), FermionicTime(key) {}
        Op(jsx::value const& jOp) : Op(jOp("flavor").int64(), jOp("key").int64()) {}
        
        jsx::value json() const {
            return jsx::object_t{
                {"key", key()},
                {"flavor", flavor()}
            };
        }
    };

    
    template<typename Mode, std::size_t type, typename Value>
    bool insert_worm_operator(Op<type> const& op, data::Data<Value> const& data, state::State<Value>& state) {
        if(!state.product().insert(op.key(), op.flavor())) return false;
        
        state.dyn().insert(op.key(), op.flavor());
           
        return true;
    };
    
    
    template<typename Mode, std::size_t type, typename Value>
    void erase_worm_operator(Op<type> const& op, data::Data<Value> const& data, state::State<Value>& state) {
        state.product().erase(op.key());
        
        state.dyn().erase(op.key(), op.flavor());
    };


    template<typename Mode, std::size_t type, typename Value>
    void init_worm_operator_data(ut::wrap<Op<type>>, jsx::value const& jParams, data::Data<Value>& data) {
    };
    
    
    using op     = Op<0>;
    using opDagg = Op<1>;
}

#endif

