#ifndef CTQMC_INCLUDE_CONFIG_TUPLE_H
#define CTQMC_INCLUDE_CONFIG_TUPLE_H


#include <tuple>

#include "worm/Op.h"
#include "worm/OpBulla.h"
#include "worm/Bilinear.h"
#include "worm/OpBullaSum.h"

#include "../Utilities.h"
#include "../../../include/JsonX.h"



namespace cfg {
    
    
    template<char const* Name, typename... Ops>
    struct WormTuple : std::tuple<Ops...> {
        static std::string name() { return Name;};  // should maybe be replaced by a cfg::worm_name<>::value ??
        
        WormTuple() = default;
        WormTuple(jsx::value const& jWormTuple) : WormTuple(jWormTuple, ut::make_sequence_t<0, sizeof...(Ops)>()) {};
        WormTuple(WormTuple const&) = default;
        WormTuple(WormTuple&&) = default;
        WormTuple& operator=(WormTuple const&) = default;
        WormTuple& operator=(WormTuple&&) = default;
        ~WormTuple() = default;
        
        template< typename... Args, typename = typename std::enable_if<!std::is_same<typename std::decay<Args>::type..., jsx::value>::value && !std::is_same<typename std::decay<Args>::type..., WormTuple>::value>::type>
        WormTuple(Args&&... args) : std::tuple<Ops...>(std::forward<Args>(args)...) {
        };
        
        template<typename... Args, typename = typename std::enable_if<!std::is_same<typename std::decay<Args>::type..., jsx::value>::value && !std::is_same<typename std::decay<Args>::type..., WormTuple>::value>::type>
        WormTuple& operator=(Args&&... args) {
            std::tuple<Ops...>::operator=(std::forward<Args>(args)...); return *this;
        };
        
        jsx::value json() const {
            return json(ut::make_sequence_t<0, sizeof...(Ops)>());
        };
        
    private:
        template<std::size_t... Indices>
        WormTuple(jsx::value const& jWorm, ut::sequence<Indices...>) : std::tuple<Ops...>(get<Indices>(jWorm)...) {};
        
        template<std::size_t... Indices>
        jsx::value json(ut::sequence<Indices...>) const {
            return jsx::object_t{
                {"name", name()},
                {"entry", jsx::array_t{std::get<Indices>(*this).json()...}}
            };
        }
        
        template<std::size_t Index>
        static jsx::value const& get(jsx::value const& jWorm) {
            return jWorm(Index);
        };
    };
    
    
    template<typename> struct worm_size;
    
    template<char const* Name, typename... Ops>
    struct worm_size<WormTuple<Name, Ops...>> {
        constexpr static std::size_t value = sizeof...(Ops);
    };
    
    template<std::size_t , typename> struct worm_op;
    
    template<std::size_t Index, char const* Name, typename... Ops>
    struct worm_op<Index, WormTuple<Name, Ops...>> {
        using type = ut::get_type_by_index_t<Index, Ops...>;
    };

    template<std::size_t Index, typename Worm> using worm_op_t = typename worm_op<Index, Worm>::type;
    
    template<typename... Args>
    std::tuple<Args...>& as_tuple(std::tuple<Args...>& arg) {
        return arg;
    }
    
    template<typename... Args>
    std::tuple<Args...> const& as_tuple(std::tuple<Args...> const& arg) {
        return arg;
    }
    
    //------------------------------------------------- insert_worm -----------------------------------------------------
    
    template<std::size_t Index, std::size_t Size, typename Worm, typename Mode, typename Value>
    struct for_each_insert {
        static bool apply(Worm const& worm, data::Data<Value> const& data, state::State<Value>& state) {
            if(!insert_worm_operator<Mode>(std::get<Index>(worm), data, state)) return false;
            return for_each_insert<Index + 1, Size, Worm, Mode, Value>::apply(worm, data, state);
        }
    };
    
    template<std::size_t Size, typename Worm, typename Mode, typename Value>
    struct for_each_insert<Size, Size, Worm, Mode, Value> {
        static bool apply(Worm const& worm, data::Data<Value> const& data, state::State<Value>& state) {
            return true;
        }
    };

    template<typename Mode, typename Value, typename... Ops>
    bool insert_worm(std::tuple<Ops...> const& worm, data::Data<Value> const& data, state::State<Value>& state) {
        return for_each_insert<0, sizeof...(Ops), std::tuple<Ops...>, Mode, Value>::apply(worm, data, state);
    };
    
    //-------------------------------------------------- erase_worm -----------------------------------------------------
    
    template<std::size_t Index, typename Worm, typename Mode, typename Value>
    struct for_each_erase {
        static void apply(Worm const& worm, data::Data<Value> const& data, state::State<Value>& state) {
            erase_worm_operator<Mode>(std::get<Index - 1>(worm), data, state);
            for_each_erase<Index - 1, Worm, Mode, Value>::apply(worm, data, state);
        };
    };
    
    template<typename Worm, typename Mode, typename Value>
    struct for_each_erase<0, Worm, Mode, Value> {
        static void apply(Worm const& worm, data::Data<Value> const& data, state::State<Value>& state) {
        };
    };
    
    template<typename Mode, typename Value, typename... Ops>
    void erase_worm(std::tuple<Ops...> const& worm, data::Data<Value> const& data, state::State<Value>& state) {
        for_each_erase<sizeof...(Ops), std::tuple<Ops...>, Mode, Value>::apply(worm, data, state);
    };
    
    //------------------------------------------------- init_worm_data --------------------------------------------------
    
    template<typename Mode, typename Value>
    struct init_worm_operator_data_functor {
        template<typename Op>
        void operator()(ut::wrap<Op> op, jsx::value const& jParams, data::Data<Value>& data) {
            init_worm_operator_data<Mode>(op, jParams, data);
        }
    };
    
    template<typename Mode, char const* Name, typename... Ops, typename Value>
    void init_worm_data(ut::wrap<WormTuple<Name, Ops...>>, jsx::value const& jParams, data::Data<Value>& data) {
        ut::for_each_type<Ops...>::apply(init_worm_operator_data_functor<Mode, Value>(), jParams, data);
    };
    
}

#endif
