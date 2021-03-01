#ifndef CTQMC_INCLUDE_CONFIG_WORM_OPBULLA_H
#define CTQMC_INCLUDE_CONFIG_WORM_OPBULLA_H


#include "Basic.h"
#include "../../impurity/Tensor.h"
#include "../../Utilities.h"
#include "../../../../include/JsonX.h"



namespace cfg {
    

    template<std::size_t type>
    struct OpBulla : Flavor, std::tuple<Flavor, Flavor, Flavor>, FermionicTime {
        OpBulla() = default;
        OpBulla(int flavor, int flavor0, int flavor1, int flavor2, ut::KeyType key = 0) : Flavor(flavor), std::tuple<Flavor, Flavor, Flavor>(flavor0, flavor1, flavor2), FermionicTime(key) {};
        OpBulla(jsx::value const& jOpBulla) : OpBulla(jOpBulla("flavor").int64(), jOpBulla("flavor0").int64(), jOpBulla("flavor1").int64(), jOpBulla("flavor2").int64(), jOpBulla("key").int64()) {};

        jsx::value json() const {
            return jsx::object_t{
                {"key", key()},
                {"flavor", flavor()},
                {"flavor0", std::get<0>(*this).flavor()},
                {"flavor1", std::get<1>(*this).flavor()},
                {"flavor2", std::get<2>(*this).flavor()}
            };
        };
    };
    

    template<typename Mode, std::size_t type, typename Value>
    bool insert_worm_operator(OpBulla<type> const& opBulla, data::Data<Value> const& data, state::State<Value>& state) {
        if(!state.product().insert(opBulla.key() - 0, std::get<0>(opBulla).flavor())) return false;
        if(!state.product().insert(opBulla.key() - 1, std::get<1>(opBulla).flavor())) return false;
        if(!state.product().insert(opBulla.key() - 2, std::get<2>(opBulla).flavor())) return false;
        
        state.dyn().insert(opBulla.key() - 0, std::get<0>(opBulla).flavor());
        state.dyn().insert(opBulla.key() - 1, std::get<1>(opBulla).flavor());
        state.dyn().insert(opBulla.key() - 2, std::get<2>(opBulla).flavor());
        
        auto const& tensor = data.template opt<imp::Tensor<Value>>();
        
        state.fact() *= tensor(opBulla.flavor()/2, std::get<0>(opBulla).flavor()/2, std::get<1>(opBulla).flavor()/2, std::get<2>(opBulla).flavor()/2);
        
        return true;
    };
    
    
    template<typename Mode, std::size_t type, typename Value>
    void erase_worm_operator(OpBulla<type> const& opBulla, data::Data<Value> const& data, state::State<Value>& state) {
        state.product().erase(opBulla.key() - 2);
        state.product().erase(opBulla.key() - 1);
        state.product().erase(opBulla.key() - 0);
        
        state.dyn().erase(opBulla.key() - 2, std::get<2>(opBulla).flavor());
        state.dyn().erase(opBulla.key() - 1, std::get<1>(opBulla).flavor());
        state.dyn().erase(opBulla.key() - 0, std::get<0>(opBulla).flavor());
        
        auto const& tensor = data.template opt<imp::Tensor<Value>>();
        
        state.fact() /= tensor(opBulla.flavor()/2, std::get<0>(opBulla).flavor()/2, std::get<1>(opBulla).flavor()/2, std::get<2>(opBulla).flavor()/2);
    };
    

    template<std::size_t type>
    bool operator<(OpBulla<type> const& lhs, OpBulla<type> const& rhs) {
        if(static_cast<Flavor const&>(lhs) < static_cast<Flavor const&>(rhs)) return true;
        if(static_cast<Flavor const&>(rhs) < static_cast<Flavor const&>(lhs)) return false;
        return static_cast<std::tuple<Flavor, Flavor, Flavor> const&>(lhs) < static_cast<std::tuple<Flavor, Flavor, Flavor> const&>(rhs);
    };

    
    template<typename Mode, typename Value>
    void init_worm_operator_data(ut::wrap<OpBulla<0>>, jsx::value const& jParams, data::Data<Value>& data) {
        auto& tensor = data.template opt<imp::Tensor<Value>>();
        if(tensor.get() == nullptr) tensor.reset(new imp::Tensor<Value>(jParams("hloc")("two body"), data.ops().flavors()/2));
    };
    
    
    using opBulla = OpBulla<0>;
}

#endif

