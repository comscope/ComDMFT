#ifndef CTQMC_INCLUDE_CONFIG_WORM_OPBULLASUM_H
#define CTQMC_INCLUDE_CONFIG_WORM_OPBULLASUM_H


#include "Basic.h"
#include "../../Utilities.h"
#include "../../impurity/Product.h"
#include "../../impurity/Observables.h"
#include "../../../../include/JsonX.h"


namespace cfg {

    
    template<std::size_t type>
    struct OpBullaSum : Flavor, FermionicTime {
        OpBullaSum() = default;
        OpBullaSum(int flavor, ut::KeyType key = 0) : Flavor(flavor), FermionicTime(key) {};
        OpBullaSum(jsx::value const& jOpBulla) : OpBullaSum(jOpBulla("flavor").int64(), jOpBulla("key").int64()) {};
        
        jsx::value json() const {
            return jsx::object_t{
                {"key", key()},
                {"flavor", flavor()}
            };
        };
    };
    
    
    template<typename Mode, std::size_t type, typename Value>
    bool insert_worm_operator(OpBullaSum<type> const& opBullaSum, data::Data<Value> const& data, state::State<Value>& state) {
        auto& product = imp::get<Mode>(state.product());
        auto const& bullaOps = imp::get<Mode>(data.template opt<imp::itf::BullaOperators<Value>>());
        
        if(product.insert(opBullaSum.key(), &bullaOps.at(opBullaSum.flavor()), opBullaSum.flavor()) == product.last()) return false;
        
        state.dyn().insert(opBullaSum.key(), opBullaSum.flavor());
        
        return true;
    };
    
    
    template<typename Mode, std::size_t type, typename Value>
    void erase_worm_operator(OpBullaSum<type> const& opBullaSum, data::Data<Value> const& data, state::State<Value>& state) {
        state.product().erase(opBullaSum.key());
        
        state.dyn().erase(opBullaSum.key(), opBullaSum.flavor());
    };

    
    template<typename Mode, std::size_t type, typename Value>
    void init_worm_operator_data(ut::wrap<OpBullaSum<type>>, jsx::value const& jParams, data::Data<Value>& data) {
        auto& bullaOps = data.template opt<imp::itf::BullaOperators<Value>>();
        if(bullaOps.get() == nullptr) bullaOps.reset(new imp::BullaOperators<Mode, Value>(jParams("hloc")("interaction"), jParams("operators"), data.eig()));
    };
    
    
    using opBullaSum = OpBullaSum<0>;
    using opBullaSumDagg = OpBullaSum<1>;
    
}

#endif

