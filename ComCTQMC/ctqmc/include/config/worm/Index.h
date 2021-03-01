
#ifndef CTQMC_INCLUDE_CONFIG_WORM_INDEX_H
#define CTQMC_INCLUDE_CONFIG_WORM_INDEX_H


#include "../Worms.h"
#include "../../Data.h"


namespace cfg {
    
    namespace worm {
    
        std::vector<int> string_to_flavors(std::vector<std::string> const& string)
        {
            std::vector<int> flavors(string.size()); auto it = flavors.begin();
            for (auto const& s : string)
                *it++ = std::stoi(s);
        
           return flavors;
        }
        
        template<typename T, typename = void> struct flavor_trait;
        
        template<typename T>
        struct flavor_trait<T, typename std::enable_if<std::is_base_of<cfg::Flavor, T>::value>::type> {
            using type = cfg::Flavor;
        };
        template<typename T>
        struct flavor_trait<T, typename std::enable_if<std::is_base_of<std::tuple<cfg::Flavor, cfg::Flavor>, T>::value>::type> {
            using type = std::tuple<cfg::Flavor, cfg::Flavor>;
        };
        
        template<typename T> using flavor_trait_t = typename flavor_trait<T>::type;

        
        template<typename Op, typename... Ops>
        struct get_dim {
            constexpr static std::size_t value = get_dim<Op>::value + get_dim<Ops...>::value;
        };
        
        template<>
        struct get_dim<cfg::Flavor> {
            constexpr static std::size_t value = 1;
        };
        
        template<>
        struct get_dim<std::tuple<cfg::Flavor, cfg::Flavor>> {
            constexpr static std::size_t value = 2;
        };

    
        template<typename> struct Index;
    
        template <>
        struct Index<cfg::partition::Worm> {
            Index(){}
            
            static std::string string(cfg::partition::Worm const& worm) {
                return "all";
            };
            
            static std::string string() {
                return "all";
            };
            
            std::size_t size() const {
                return 1;
            };
            
            std::size_t integer(cfg::partition::Worm const& worm) const {
                return 0;
            };

        };
    
        
        
        template<char const* Name, typename... Ops>
        struct Index<cfg::WormTuple<Name, Ops...>> {
            Index() = delete;
            Index(std::size_t const N) :
            N_(N), D_(std::pow(N_, get_dim<flavor_trait_t<Ops>...>::value)) {
            };
            
            static std::string string(cfg::WormTuple<Name, Ops...> const& worm) {
                return get_string(worm, ut::make_sequence_t<0, sizeof...(Ops)>());
            };
            
            std::size_t size() const {
                return D_;
            };
            
            std::size_t integer(cfg::WormTuple<Name, Ops...> const& worm) const {
                return get_integer(worm, ut::make_sequence_t<0, sizeof...(Ops)>());
            };
            
        private:
            std::size_t const N_, D_;
            
            template<std::size_t... Indices>
            static std::string get_string(cfg::WormTuple<Name, Ops...> const& worm, ut::sequence<Indices...>) {
                return string_impl(std::get<Indices>(worm)...);
            }
            
            static std::string string_impl(cfg::Flavor const& op) {
                return std::to_string(op.flavor());
            };
            
            template<typename... Args>
            static std::string string_impl(std::tuple<cfg::Flavor, cfg::Flavor> const& ops, Args const&... args) {
                return string_impl(std::get<0>(ops)) + "_" + string_impl(std::get<1>(ops), args...);
            };
            
            template<typename... Args>
            static std::string string_impl(cfg::Flavor const& op, Args const&... args) {
                return  string_impl(op) + "_" + string_impl(args...);
            };
            
            //Get integer from worm
            template<std::size_t... Indices>
            std::size_t get_integer(cfg::WormTuple<Name, Ops...> const& worm, ut::sequence<Indices...>) const {
                return integer_impl(std::get<Indices>(worm)...);
            };
            
            std::size_t integer_impl(cfg::Flavor const& op) const {
                return op.flavor();
            };
            
            template<typename... Args>
            std::size_t integer_impl(std::tuple<cfg::Flavor, cfg::Flavor> const& ops, Args const&... args) const {
                return integer_impl(std::get<0>(ops)) + N_*integer_impl(std::get<1>(ops), args...);
            };
            
            template<typename... Args>
            std::size_t integer_impl(cfg::Flavor const& op, Args const&... args) const {
                return  integer_impl(op) + N_*integer_impl(args...);
            };
            
        };
        
    }
    
}

#endif
