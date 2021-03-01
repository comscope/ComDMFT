#ifndef CTQMC_INCLUDE_CONFIG_VARIANT_H
#define CTQMC_INCLUDE_CONFIG_VARIANT_H

#include <stdexcept>
#include <iostream>
#include <array>
#include <vector>
#include <algorithm>

#include "../Utilities.h"
#include "../../../include/JsonX.h"

namespace cfg {
    //----------------------------------------------------------------------------------------------------------------------------
    
    // Dae non-trivial-constructor-destructor terror muess besser tested werde !!   Allgmein, da gits glaub no einiges z'verbessere ......

    namespace impl {
        
        template<typename... Args> union variant_union {};
        
        template<typename T, typename... Args>
        union variant_union<T, Args...> {
            T value;
            variant_union<Args...> next;
            variant_union() {};
            ~variant_union() {};
        };
        
        
        template<typename T, typename U, typename... Args>
        struct get_by_type {
            static T& at(variant_union<U, Args...>& variant) {
                return get_by_type<T, Args...>::at(variant.next);
            };
            static T const& at(variant_union<U, Args...> const& variant) {
                return get_by_type<T, Args...>::at(variant.next);
            };
        };
        
        template<typename T, typename... Args>
        struct get_by_type<T, T, Args...> {
            static T& at(variant_union<T, Args...>& variant) {
                return variant.value;
            };
            static T const& at(variant_union<T, Args...> const& variant) {
                return variant.value;
            };
        };

        
        template<typename... Args> struct apply_by_index;
        
        template<typename T, typename... Args>
        struct apply_by_index<T, Args...> {
            template<typename V, typename F, typename... FArgs>
            static void to(int index, V&& variant, F&& func, FArgs&&... fargs) {
                if(index == 0)
                    func(variant.value, std::forward<FArgs>(fargs)...);
                else
                    apply_by_index<Args...>::to(index - 1, variant.next, std::forward<F>(func), std::forward<FArgs>(fargs)...);
            }
        };
        
        template<>
        struct apply_by_index<> {
            template<typename V, typename F, typename... FArgs>
            static void to(int index, V&& variant, F&& func, FArgs&&... fargs) {
                throw std::runtime_error("cfg::impl::apply_by_index: invalid index");
            }
        };

        
        struct construct {
            template<typename T, typename... Args>
            void operator()(T& t, Args&&... args) const {
                new(&t) T(std::forward<Args>(args)...);
            }
        };
        
        struct write {
            template<typename T>
            void operator()(T const& t, jsx::value& jWorm) const {
                jWorm = t.json();
            }
        };
        
        struct destruct {
            template<typename T, typename... Args>
            void operator()(T& t, Args&&... args) const {
                t.~T();
            }
        };

    }

    
    template<typename... Args>
    struct WormVariant {
        constexpr static int size() {
            return sizeof...(Args);
        };
        
        std::vector<std::string> static get_names() {
            return std::vector<std::string>{Args::name()...};
        };
        
        WormVariant() = delete;
        WormVariant(jsx::value jWorm) {
            auto names = get_names();
            
            auto it = std::find(names.begin(), names.end(), jWorm("name").string());
            if(it == names.end())
                throw std::runtime_error("cfg::WormVariant::WormVariant: worm with name " + jWorm("name").string() + " not found !");

            apply_by_index::to(index_ = it - names.begin(), variant_, impl::construct(), jWorm("entry"));
        }
        WormVariant(WormVariant const&) = delete;
        WormVariant(WormVariant&&) = delete;
        WormVariant& operator=(WormVariant const&) = delete;
        WormVariant& operator=(WormVariant&&) = delete;
        ~WormVariant() {
            apply_by_index::to(index_, variant_, impl::destruct());
        };
        
        template<typename T, typename std::enable_if<!std::is_same<typename std::decay<T>::type, WormVariant>::value, int>::type = 0>
        T& operator=(T&& arg) {
            using Entry = typename std::decay<T>::type;
            
            apply_by_index::to(index_, variant_, impl::destruct());
            index_ = get_index<Entry>::value;
            new(&get_by_type<Entry>::at(variant_)) Entry(std::forward<T>(arg));

            return get_by_type<Entry>::at(variant_);
        }
        
        std::size_t index() const {
            return index_;
        }
        
        jsx::value json() const {
            jsx::value jWorm;
            apply_by_index::to(index_, variant_, impl::write(), jWorm);
            return jWorm;
        }
        
    private:
        template<typename T> using get_by_type = impl::get_by_type<T, Args...>;
        template<typename T> using get_index = ut::get_index_by_type<T, Args...>;
        using apply_by_index = impl::apply_by_index<Args...>;
        
        int index_; 
        impl::variant_union<Args...> variant_;
        
       
        template<typename T, typename... OtherArgs> friend T& get(WormVariant<OtherArgs...>&);
        template<typename T, typename... OtherArgs> friend T const& get(WormVariant<OtherArgs...> const&);
        template<typename... OtherArgs, typename F, typename... FArgs> friend void apply(WormVariant<OtherArgs...>&, F&&, FArgs&&...);
        template<typename... OtherArgs, typename F, typename... FArgs> friend void apply(WormVariant<OtherArgs...> const&, F&&, FArgs&&...);
        
    };
    
    
    template<typename... Args> struct get_index;
    
    template<typename T, typename... Args>
    struct get_index<T, WormVariant<Args...>> {
        static constexpr int value = ut::get_index_by_type<T, Args...>::value;
    };

    
    template<typename... Args> struct for_each_type;
    
    template<typename... Args>
    struct for_each_type<WormVariant<Args...>> {
        template<typename F, typename... FArgs>
        static void apply(F&& func, FArgs&&... fargs) {
            ut::for_each_type<Args...>::apply(std::forward<F>(func), std::forward<FArgs>(fargs)...);
        }
    };
    
    
    template<typename T, typename... Args>
    T& get(WormVariant<Args...>& variant) {
        if(ut::get_index_by_type<T, Args...>::value != variant.index_)
            throw std::runtime_error("cfg::get: access of inactive variant");
        return impl::get_by_type<T, Args...>::at(variant.variant_);
    }
    
    template<typename T, typename... Args>
    T const& get(WormVariant<Args...> const& variant) {
        if(ut::get_index_by_type<T, Args...>::value != variant.index_)
            throw std::runtime_error("cfg::get: access of inactive variant");
        return impl::get_by_type<T, Args...>::at(variant.variant_);
    }
    
    
    template<typename... Args, typename F, typename... FArgs>
    void apply(WormVariant<Args...>& variant, F&& func, FArgs&&... fargs) {
        impl::apply_by_index<Args...>::to(variant.index_, variant.variant_, std::forward<F>(func), std::forward<FArgs>(fargs)...);
    }
    
    template<typename... Args, typename F, typename... FArgs>
    void apply(WormVariant<Args...> const& variant, F&& func, FArgs&&... fargs) {
        impl::apply_by_index<Args...>::to(variant.index_, variant.variant_, std::forward<F>(func), std::forward<FArgs>(fargs)...);
    }
    
}

#endif
