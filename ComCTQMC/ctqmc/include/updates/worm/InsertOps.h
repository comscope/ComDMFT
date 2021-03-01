#ifndef CTQMC_INCLUDE_UPDATES_WORM_INSERTOPS_H
#define CTQMC_INCLUDE_UPDATES_WORM_INSERTOPS_H

#include <tuple>

#include "Ops.h"
#include "../../Data.h"
#include "../../State.h"


namespace upd {
    
    namespace worm {
        
        
        template<typename... Ops>
        struct ins_functor {
            ins_functor(std::vector<double>& prob, std::vector<std::tuple<Ops...>>& candidates) : distr_(prob), candidates_(candidates) {};
            ~ins_functor() {
                for(auto& p : distr_) p /= distr_.back();
            };
            
            void operator()(std::tuple<Ops...> candidate) {
                distr_.push_back((distr_.size() ? distr_.back() : .0) + 1.);
                candidates_.push_back(candidate);
            }
            
        private:
            std::vector<double>& distr_;
            std::vector<std::tuple<Ops...>>& candidates_;
        };
        
        
        template<std::size_t Counter>
        struct set_keys {
            template<typename... Ops>
            static void apply(std::tuple<Ops...>& candidate, ut::UniformRng& urng) {
                std::get<Counter - 1>(candidate).key() = ut::KeyMax*urng();
                set_keys<Counter - 1>::apply(candidate, urng);
            }
        };
        
        template<>
        struct set_keys<0> {
            template<typename... Ops>
            static void apply(std::tuple<Ops...>& candidate, ut::UniformRng& urng) {
                std::get<0>(candidate).key() = ut::KeyMax*urng();
            }
        };
        
        
        template<typename...> struct InsertOps;
        
        
        template<typename Origin, typename... Ops, typename Target, std::size_t... Indices>
        struct InsertOps<Origin, std::tuple<Ops...>, Target, ut::sequence<Indices...>> {
            template<typename Value>
            InsertOps(jsx::value const& jOpsList, data::Data<Value> const& data) {
                ins_functor<Ops...> functor(distr_, candidates_);
                
                if(jOpsList.is<jsx::empty_t>())
                    all_ops<Ops...>::apply(data, functor);
                else
                    for(auto const& jOp : jOpsList.array())
                        read_op<Ops...>::apply(data, jOp.array().begin(), functor);
            }
            
            template<typename Value>
            bool propose(double const urn, data::Data<Value> const& data, state::State<Value>& state, ut::UniformRng& urng) {
                auto it = std::upper_bound(distr_.begin(), distr_.end(), urn*distr_.back());
                prob_ = (it != distr_.begin() ? *it - *(it - 1) : *it);
                insert_ = candidates_[it - distr_.begin()];
                
                set_keys<sizeof...(Ops)>::apply(insert_, urng);

                target_ = std::make_tuple(std::get<Indices>(std::tuple_cat(cfg::as_tuple(cfg::get<Origin>(state.worm())), insert_))...);
                return true;
            };
            
            Target const& target(Origin const& origin) const {
                return target_;
            };
            
            template<typename Mode, typename Value>
            bool impurity(data::Data<Value> const& data, state::State<Value>& state) const {
                return cfg::insert_worm<Mode>(insert_, data, state);
            };
            
            template<typename Value>
            double ratio(data::Data<Value> const& data, state::State<Value>& state) {
                return std::pow(ut::beta(), sizeof...(Ops))/prob_;
            };
            
            template<typename Value>
            void accept(data::Data<Value> const& data, state::State<Value>& state) {
                state.worm() = target_;
            };
            
            template<typename Value>
            void reject(data::Data<Value> const& data, state::State<Value>& state) {
            };
            
            jsx::value json() const {
                return jsx::null_t(); // TODO
            };
            
        private:
            
            std::vector<double> distr_;
            std::vector<std::tuple<Ops...>> candidates_;
            
            double prob_;
            
            Target target_;
            std::tuple<Ops...> insert_;

        };


    }
}

#endif
