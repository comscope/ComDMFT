#ifndef CTQMC_INCLUDE_UPDATES_WORM_OPS_H
#define CTQMC_INCLUDE_UPDATES_WORM_OPS_H

#include "../include/Surviving.h"
#include "../../config/Worms.h"

namespace upd {
    
    namespace worm {
        
        
        template<typename...> struct all_ops;
        
        template<std::size_t type, typename... Ops>
        struct all_ops<cfg::Op<type>, Ops...> {
            template<typename Value, typename... FArgs>
            static void apply(data::Data<Value> const& data, FArgs&&... fargs) {
                for(int f = type; f < data.ops().flavors(); f += 2)
                    all_ops<Ops...>::apply(data, std::forward<FArgs>(fargs)..., cfg::Op<type>(f));
            }
        };
        
        template<std::size_t type1, std::size_t type2,  typename... Ops>
        struct all_ops<cfg::Bilinear<type1, type2>, Ops...> {
            template<typename Value, typename... FArgs>
            static void apply(data::Data<Value> const& data, FArgs&&... fargs) {
                for(int f1 = type1; f1 < data.ops().flavors(); f1 += 2)
                    for(int f2 = type2; f2 < data.ops().flavors(); f2 += 2)
                        if(f1 != f2) all_ops<Ops...>::apply(data, std::forward<FArgs>(fargs)..., cfg::Bilinear<type1, type2>(f1, f2));
            }
        };
        
        template<std::size_t type, typename... Ops>
        struct all_ops<cfg::OpBulla<type>, Ops...> {
            template<typename Value, typename... FArgs>
            static void apply(data::Data<Value> const& data, FArgs&&... fargs) {
                for(int f = type; f < data.ops().flavors(); f += 2) {
                    auto const& tensor = data.template opt<imp::Tensor<Value>>();
                    for(auto const& nonzero : tensor.nonzero(f/2))
                        all_ops<Ops...>::apply(data, std::forward<FArgs>(fargs)..., cfg::OpBulla<type>{f, nonzero.flavor1(), nonzero.flavor2(), nonzero.flavor3()});
                }
            }
        };
        
        template<std::size_t type, typename... Ops>
        struct all_ops<cfg::OpBullaSum<type>, Ops...> {
            template<typename Value, typename... FArgs>
            static void apply(data::Data<Value> const& data, FArgs&&... fargs) {
                for(int f = type; f < data.ops().flavors(); f += 2)
                    all_ops<Ops...>::apply(data, std::forward<FArgs>(fargs)..., cfg::OpBullaSum<type>(f));
            }
        };
        
        template<>
        struct all_ops<> {
            template<typename Value, typename F, typename... Args>
            static void apply(data::Data<Value> const& data, F&& func, Args&&... args) {
                if(surviving(data, std::forward<Args>(args)...))
                    func(std::make_tuple(std::forward<Args>(args)...));
            }
        };
        
        
        
        
        template<typename...> struct read_op;
        
        template<std::size_t type, typename... Ops>
        struct read_op<cfg::Op<type>, Ops...> {
            template<typename Value, typename... FArgs>
            static void apply(data::Data<Value> const& data, jsx::array_t::const_iterator it, FArgs&&... fargs) {
                read_op<Ops...>::apply(data, it + 1, std::forward<FArgs>(fargs)..., cfg::Op<type>(it->int64()));
            };
        };
        
        template<std::size_t type1, std::size_t type2, typename... Ops>
        struct read_op<cfg::Bilinear<type1, type2>, Ops...> {
            template<typename Value, typename... FArgs>
            static void apply(data::Data<Value> const& data, jsx::array_t::const_iterator it, FArgs&&... fargs) {
                read_op<Ops...>::apply(data, it + 2, std::forward<FArgs>(fargs)..., cfg::Bilinear<type1, type2>(it->int64(), (it + 1)->int64()));
            };
        };
        
        template<std::size_t type, typename... Ops>
        struct read_op<cfg::OpBulla<type>, Ops...> {
            template<typename Value, typename... FArgs>
            static void apply(data::Data<Value> const& data, jsx::array_t::const_iterator it, FArgs&&... fargs) {
                auto const& tensor = data.template opt<imp::Tensor<Value>>();
                for(auto const& nonzero : tensor.nonzero(it->int64()/2))
                    read_op<Ops...>::apply(data, it + 1, std::forward<FArgs>(fargs)..., cfg::OpBulla<type>(it->int64(), nonzero.flavor1(), nonzero.flavor2(), nonzero.flavor3()));
            };
        };
        
        template<std::size_t type, typename... Ops>
        struct read_op<cfg::OpBullaSum<type>, Ops...> {
            template<typename Value, typename... FArgs>
            static void apply(data::Data<Value> const& data, jsx::array_t::const_iterator it, FArgs&&... fargs) {
                read_op<Ops...>::apply(data, it + 1, std::forward<FArgs>(fargs)..., cfg::OpBullaSum<type>(it->int64()));
            };
        };
        
        template<>
        struct read_op<> {
            template<typename Value, typename F, typename... Args>
            static void apply(data::Data<Value> const& data, jsx::array_t::const_iterator it, F&& func, Args&&... args) {
                func(std::make_tuple(std::forward<Args>(args)...));
            };
        };
        
        
    }
    
}

#endif
