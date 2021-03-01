#ifndef CTQMC_INCLUDE_UPDATES_WORM_RECONNECT_H
#define CTQMC_INCLUDE_UPDATES_WORM_RECONNECT_H

#include <tuple>

#include "../include/Generic.h"
#include "../../bath/algebra/Exchange.h"
#include "../../Data.h"
#include "../../State.h"

#include <iostream>

namespace upd {
    
    namespace worm {
        
        template<typename Worm, std::size_t Index, typename = cfg::worm_op_t<Index, Worm>> struct Reconnect;
        
        
        //---------------------------------------------------------- Op ----------------------------------------------------------------

        template<char const* name, typename... Ops, std::size_t Index, std::size_t type>
        struct Reconnect<cfg::WormTuple<name, Ops...>, Index, cfg::Op<type>> {
            using Op   = cfg::Op<type>;
            using Worm = cfg::WormTuple<name, Ops...>;
            
            Reconnect() = default;
            
            template<typename Value>
            bool propose(double const urn, data::Data<Value> const& data, state::State<Value>& state, ut::UniformRng& urng) {
                target_ = cfg::get<Worm>(state.worm());

                auto& ops = state.expansion()[std::get<Index>(target_).flavor()];
                if(!ops.size()) return false;
                
                connect_ = std::get<Index>(target_);
                std::get<Index>(target_).key() = ops[urng()*ops.size()];
                
                return true;
            };
            
            Worm const& target(Worm const& origin) const {
                return target_;
            };
            
            template<typename Value>
            double ratio(data::Data<Value> const& data, state::State<Value>& state) {
                return 1.;
            };
            
            template<typename Value>
            void bath(data::Data<Value> const& data, std::vector<bath::Bath<Value>>& baths) const {
                baths[data.hyb().bath(connect_.flavor())].add(data.hyb(), bath::Exchange{std::get<Index>(target_).key(), connect_.key(), connect_.flavor()});
            };
            
            template<typename Value>
            void accept(data::Data<Value> const& data, state::State<Value>& state) {
                auto& ops = state.expansion()[connect_.flavor()];
                ops.erase(std::get<Index>(target_).key()); ops.insert(connect_.key());
                
                state.signTimeOrder() *= -1;
                
                state.worm() = target_;
            };
            
            template<typename Value>
            void reject(data::Data<Value> const& data, state::State<Value>& state) {
            };
            
            jsx::value json() const {
                return jsx::null_t(); // TODO
            };
            
        private:
            Worm target_;
            Op connect_;
        };
        
        
        //---------------------------------------------------------- Bilinear ----------------------------------------------------------------
        
        template<char const* name, typename... Ops, std::size_t Index, std::size_t... types>
        struct Reconnect<cfg::WormTuple<name, Ops...>, Index, cfg::Bilinear<types...>> {
            using Bilinear = cfg::Bilinear<types...>;
            using Worm     = cfg::WormTuple<name, Ops...>;
            
            Reconnect() = default;
            
            template<typename Value>
            bool propose(double const urn, data::Data<Value> const& data, state::State<Value>& state, ut::UniformRng& urng) {
                target_ = cfg::get<Worm>(state.worm());  new_ = old_ = std::get<Index>(target_);
                
                flavor0_ = std::get<0>(old_).flavor();  auto const& ops0 = state.expansion()[flavor0_];
                flavor1_ = std::get<1>(old_).flavor();  auto const& ops1 = state.expansion()[flavor1_];
                
                if(!(ops0.size() && ops1.size())) return false;
                
                auto join = upper_bound(ops0, flavor0_, ops1, flavor1_, new_.key() = ops1[urng()*ops1.size()]);
                
                if(join.second != flavor0_) return false;  
                
                new_.key() += 1;
                
                auto splitHigh = upper_bound(ops0, flavor0_, ops1, flavor1_, old_.key());
                
                if(join.first == splitHigh.first) return false;
                
                auto joinHigh = upper_bound(ops0, flavor0_, ops1, flavor1_, join.first);
                
                if(ut::cyclic(joinHigh.first - join.first) > ut::cyclic(old_.key() - 1 - join.first)) joinHigh.first = old_.key() - 1;
                
                diffSplit_ = ut::cyclic(splitHigh.first - old_.key());
                diffJoin_  = ut::cyclic(joinHigh.first  - new_.key());
                
                keySplit_ = ut::cyclic(splitHigh.first - 1 - urng()*diffSplit_);
                keyJoin_  = join.first;  std::get<Index>(target_) = new_;
                
                return true;
            };
            
            Worm const& target(Worm const& origin) const {
                return target_;
            };
            
            template<typename Value>
            double ratio(data::Data<Value> const& data, state::State<Value>& state) {
                return static_cast<double>(diffSplit_)/static_cast<double>(diffJoin_);
            };
            
            template<typename Value>
            void bath(data::Data<Value> const& data, std::vector<bath::Bath<Value>>& baths) const {
                baths[data.hyb().bath(flavor0_)].add(data.hyb(), bath::Exchange{      keyJoin_,      keySplit_, flavor0_});
                baths[data.hyb().bath(flavor1_)].add(data.hyb(), bath::Exchange{new_.key() - 1, old_.key() - 1, flavor1_});
            };
            
            
            template<typename Mode, typename Value>
            bool impurity(data::Data<Value> const& data, state::State<Value>& state) const {
                state.product().erase(old_.key());
                if(!state.product().insert(new_.key(), flavor0_)) return false;
                
                state.product().erase(keyJoin_);
                if(!state.product().insert(keySplit_, flavor0_)) return false;
                
                state.dyn().erase(old_.key(), flavor0_);
                state.dyn().insert(new_.key(), flavor0_);
                
                state.dyn().erase(keyJoin_, flavor0_);
                state.dyn().insert(keySplit_, flavor0_);
                
                return true;
            };
            
            
            template<typename Value>
            void accept(data::Data<Value> const& data, state::State<Value>& state) {
                state.expansion()[flavor0_].erase(keyJoin_);
                state.expansion()[flavor1_].erase(new_.key() - 1);
                
                state.expansion()[flavor0_].insert(keySplit_);
                state.expansion()[flavor1_].insert(old_.key() - 1);
                
                state.signTimeOrder() *= -1;
                
                state.worm() = target_;
            };
            
            template<typename Value>
            void reject(data::Data<Value> const& data, state::State<Value>& state) {
            };
            
            jsx::value json() const {
                return jsx::null_t(); // TODO
            };
            
        private:
            
            ut::KeyType keyJoin_, keySplit_;
            ut::KeyType diffJoin_, diffSplit_;
            int flavor0_, flavor1_;
            Bilinear new_, old_;
            Worm target_;
            
            std::pair<ut::KeyType, int> upper_bound(cfg::Entries const& ops0, int flavor0, cfg::Entries const& ops1, int flavor1, ut::KeyType key) {
                auto it0 = std::upper_bound(ops0.begin(), ops0.end(), key);
                auto it1 = std::upper_bound(ops1.begin(), ops1.end(), key);
                
                ut::KeyType key0 = it0 != ops0.end() ? *it0 : *(it0 = ops0.begin()) + ut::KeyMax;
                ut::KeyType key1 = it1 != ops1.end() ? *it1 : *(it1 = ops1.begin()) + ut::KeyMax;
                
                return key0 > key1 ? std::make_pair(*it1, flavor1) : std::make_pair(*it0, flavor0);
            }
        };
        
    }
    
}

#endif
