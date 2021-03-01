#ifndef CTQMC_INCLUDE_OBSERVABLES_WORM_MEASUREMENT_H
#define CTQMC_INCLUDE_OBSERVABLES_WORM_MEASUREMENT_H

#include "Basis.h"
#include "../../Data.h"


namespace obs {
    
    namespace worm {
        
        enum class MeasType { Static, Dynamic };
    
        
        template<typename T, typename = void> struct time_trait;
        
        template<typename T>
        struct time_trait<T, typename std::enable_if<std::is_base_of<cfg::FermionicTime, T>::value>::type> {
            using type = cfg::FermionicTime;
        };
        template<typename T>
        struct time_trait<T, typename std::enable_if<std::is_base_of<cfg::BosonicTime, T>::value>::type> {
            using type = cfg::BosonicTime;
        };
        
        template<typename T> using time_trait_t = typename time_trait<T>::type;
        
        
        
        template<typename, typename, FuncType, MeasType, typename> struct Meas;
        
        template<typename Mode, typename Value, FuncType funcType, char const* Name, typename... Ops>
        struct Meas<Mode, Value, funcType, MeasType::Static, cfg::WormTuple<Name, Ops...>> : Basis<Value, funcType, time_trait_t<Ops>...> {
            Meas() = delete;
            Meas(jsx::value const& jWorm, data::Data<Value> const& data) : basis_type(jWorm, data) {
            };
            Meas(Meas const&) = delete;
            Meas(Meas&&) = default;
            Meas& operator=(Meas const&) = delete;
            Meas& operator=(Meas&&) = delete;
            ~Meas() = default;
            
            void add(Value const sign, data::Data<Value> const& data, state::State<Value>& state) {
                add_impl(sign, cfg::get<worm_type>(state.worm()), ut::make_sequence_t<0, sizeof...(Ops)>());
            }
            
        private:
            using worm_type = cfg::WormTuple<Name, Ops...>;
            using basis_type = Basis<Value, funcType, time_trait_t<Ops>...>;
            
            template<std::size_t... Indices>
            void add_impl(Value const sign, worm_type const& worm, ut::sequence<Indices...>) {
                basis_type::add(sign, std::get<Indices>(worm)...);
            }
        };
        
        
        template<typename Mode, typename Value, FuncType funcType, char const* Name, std::size_t opType, typename... Ops>
        struct Meas<Mode, Value, funcType, MeasType::Dynamic, cfg::WormTuple<Name, cfg::Op<opType>, Ops...>> : Basis<Value, funcType, time_trait_t<cfg::Op<opType>>, time_trait_t<Ops>...> {
            Meas() = delete;
            Meas(jsx::value const& jWorm, data::Data<Value> const& data) : basis_type(jWorm, data) {
            };
            Meas(Meas const&) = delete;
            Meas(Meas&&) = default;
            Meas& operator=(Meas const&) = delete;
            Meas& operator=(Meas&&) = delete;
            ~Meas() = default;
            
            void add(Value const sign, data::Data<Value> const& data, state::State<Value>& state) {
                add_impl(sign, *data.dyn(), state.dyn(), imp::get<Mode>(state.densityMatrix()), cfg::get<worm_type>(state.worm()), ut::make_sequence_t<0, 1 + sizeof...(Ops)>());
            }
            
        private:
            using worm_type = cfg::WormTuple<Name, cfg::Op<opType>, Ops...>;
            using basis_type = Basis<Value, funcType, time_trait_t<cfg::Op<opType>>, time_trait_t<Ops>...>;
            
            template<std::size_t... Indices>
            void add_impl(Value const sign, imp::Simple const& dynFunc, imp::itf::Dynamic const& dynWeight, imp::DensityMatrix<Mode, Value> const& densityMatrix, worm_type const& worm, ut::sequence<Indices...>) {
                Value boundary = .0;
                
                for(auto sec : densityMatrix)
                    boundary += densityMatrix.weight(sec)*dynFunc.D0Qf(sec, std::get<0>(worm).flavor());
                
                basis_type::add(sign*(boundary - dynWeight.fkinks(std::get<0>(worm).flavor(), std::get<0>(worm).key())), std::get<Indices>(worm)...);
            }
        };
        
        
        template<typename Mode, typename Value, FuncType funcType, char const* Name, typename... Ops>
        struct Meas<Mode, Value, funcType, MeasType::Dynamic, cfg::WormTuple<Name, Ops...>> {
            Meas() = delete;
            Meas(jsx::value const& jWorm, data::Data<Value> const& data) {
                throw std::runtime_error("obs::worm::Meas: MeasType::Dynamic for " + std::string(Name) + " worm makes no sense");
            };
            Meas(Meas const&) = delete;
            Meas(Meas&&) = default;
            Meas& operator=(Meas const&) = delete;
            Meas& operator=(Meas&&) = delete;
            ~Meas() = default;
            
            void add(Value const sign, data::Data<Value> const& data, state::State<Value>& state) {
            };
            
            void store(jsx::value& measurements, std::int64_t samples) {
            };
            
        };
    }
    
}

#endif
