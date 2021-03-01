#ifndef CTQMC_INCLUDE_UPDATES_WORM_REMOVEOPS_H
#define CTQMC_INCLUDE_UPDATES_WORM_REMOVEOPS_H

#include <tuple>

#include "Ops.h"
#include "../../Data.h"
#include "../../State.h"


namespace upd {
    
    namespace worm {
        
        template<typename... Ops>
        struct rem_functor {
            rem_functor(std::map<std::tuple<Ops...>, double>& candidateprobs) : norm_(.0), candidateprobs_(candidateprobs) {};
            ~rem_functor() {
                for(auto& x : candidateprobs_) x.second /= norm_;
            };
            
            void operator()(std::tuple<Ops...> candidate) {
                if(candidateprobs_.find(candidate) != candidateprobs_.end())
                    throw std::runtime_error("upd::worm::rem_str: candidate appears twice");
                norm_ += 1.; candidateprobs_[candidate] = 1.;
            }
            
        private:
            double norm_;
            std::map<std::tuple<Ops...>, double>& candidateprobs_;
        };
        
        
        template<typename...> struct RemoveOps;
        

        template<typename Origin, typename Target, typename... Ops, std::size_t... Indices>
        struct RemoveOps<Origin, Target, std::tuple<Ops...>, ut::sequence<Indices...>> {
            using target_sequence = ut::transform_sequence_t<ut::sequence<Indices...>, ut::make_sequence_t<0, cfg::worm_size<Target>::value>>;
            using remove_sequence = ut::transform_sequence_t<ut::sequence<Indices...>, ut::make_sequence_t<cfg::worm_size<Target>::value, cfg::worm_size<Origin>::value>>;
            
            template<typename Value>
            RemoveOps(jsx::value const& jOpsList, data::Data<Value> const& data) {
                rem_functor<Ops...> functor(candidateprobs_);
                
                if(jOpsList.is<jsx::empty_t>())
                    all_ops<Ops...>::apply(data, functor);
                else
                    for(auto const& jOp : jOpsList.array())
                        read_op<Ops...>::apply(data, jOp.array().begin(), functor);

            }

            template<typename Value>
            bool propose(double urn, data::Data<Value> const& data, state::State<Value>& state, ut::UniformRng& urng) {
                remove_ = select(cfg::get<Origin>(state.worm()), remove_sequence());
                
                auto it = candidateprobs_.find(remove_);
                if(it != candidateprobs_.end())
                    prob_ = it->second;
                else
                    return false;
                
                target_ = select(cfg::get<Origin>(state.worm()), target_sequence());
                
                return true;
            };
            
            Target const& target(Origin const& origin) const {
                return target_;
            };
            
            template<typename Mode, typename Value>
            bool impurity(data::Data<Value> const& data, state::State<Value>& state) const {
                cfg::erase_worm<Mode>(remove_, data, state);  return true;
            };
            
            template<typename Value>
            double ratio(data::Data<Value> const& data, state::State<Value>& state) {
                return prob_/std::pow(ut::beta(), sizeof...(Ops));
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
            
        protected:
            Target target_;
            std::tuple<Ops...> remove_;
            
            std::map<std::tuple<Ops...>, double> candidateprobs_;
            double prob_;

            template<std::size_t... Select>
            static auto select(Origin const& origin, ut::sequence<Select...>)
            -> decltype(std::make_tuple(std::get<Select>(cfg::as_tuple(origin))...)) {
                return std::make_tuple(std::get<Select>(cfg::as_tuple(origin))...);
            }
        };

    }
}

#endif
