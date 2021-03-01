#ifndef CTQMC_INCLUDE_CONFIG_WORM_BILINEAR_H
#define CTQMC_INCLUDE_CONFIG_WORM_BILINEAR_H

#include "Basic.h"
#include "../../Utilities.h"
#include "../../../../include/JsonX.h"


namespace cfg {
    
    template<std::size_t type1, std::size_t type2>
    struct Bilinear : std::tuple<Flavor, Flavor>, BosonicTime {
        Bilinear() = default;
        Bilinear(int flavor0, int flavor1, ut::KeyType key = 0) : std::tuple<Flavor, Flavor>(flavor0, flavor1), BosonicTime(key) {}
        Bilinear(jsx::value const& jBilinear) : Bilinear(jBilinear("flavor0").int64(), jBilinear("flavor1").int64(), jBilinear("key").int64()) {}
        
        jsx::value json() const {
            return jsx::object_t{
                {"key", key()},
                {"flavor0", std::get<0>(*this).flavor()},
                {"flavor1", std::get<1>(*this).flavor()}
            };
        }
    };
    
    
    template<typename Mode, std::size_t... types, typename Value>
    bool insert_worm_operator(Bilinear<types...> const& bilinear, data::Data<Value> const& data, state::State<Value>& state) {
        if(!state.product().insert(bilinear.key() - 0, std::get<0>(bilinear).flavor())) return false;
        if(!state.product().insert(bilinear.key() - 1, std::get<1>(bilinear).flavor())) return false;
        
        state.dyn().insert(bilinear.key() - 0, std::get<0>(bilinear).flavor());
        state.dyn().insert(bilinear.key() - 1, std::get<1>(bilinear).flavor());
        
        return true;
    };
    
    
    template<typename Mode, std::size_t... types, typename Value>
    void erase_worm_operator(Bilinear<types...> const& bilinear, data::Data<Value> const& data, state::State<Value>& state) {
        state.product().erase(bilinear.key() - 1);
        state.product().erase(bilinear.key() - 0);
        
        state.dyn().erase(bilinear.key() - 1, std::get<1>(bilinear).flavor());
        state.dyn().erase(bilinear.key() - 0, std::get<0>(bilinear).flavor());
    };

    
    template<typename Mode, std::size_t... types, typename Value>
    void init_worm_operator_data(ut::wrap<Bilinear<types...>>, jsx::value const& jParams, data::Data<Value>& data) {
    };
    
    
    using bilinearPH = Bilinear<1, 0>;
    using bilinearPP = Bilinear<1, 1>;
    using bilinearHH = Bilinear<0, 0>;

}

#endif

