#ifndef CTQMC_INCLUDE_UPDATES_INCLUDE_GENERIC_H
#define CTQMC_INCLUDE_UPDATES_INCLUDE_GENERIC_H


#include "ImplementBath.h"
#include "ImplementImpurity.h"

#include "../../markovchain/Update.h"

#include "../../Data.h"
#include "../../State.h"


namespace upd {
    
    template<typename...> struct get_target_func;
    
    template<typename T, typename U, typename O>
    struct get_target_func<T const& (U::*)(O const&) const> {
        using Target = T;
        using Origin = O;
    };
    
    
    template<typename Upd, typename Mode, typename Value>
    struct Generic : mch::itf::Update<Value> {
        using Origin = typename get_target_func<decltype(&Upd::target)>::Origin;
        using Target = typename get_target_func<decltype(&Upd::target)>::Target;
        
        template<typename Space> using get_index = cfg::get_index<Space, cfg::Worm>;
        
        Generic() = delete;
        template<typename... Args>
        Generic(double probChoose, Args&&... args) : mch::itf::Update<Value>(get_index<Origin>::value, get_index<Target>::value, probChoose),
        flag_(ut::Flag::Reject), update_(std::forward<Args>(args)...) {
        };
        Generic(Generic const&) = delete;
        Generic(Generic&&) = delete;
        Generic& operator=(Generic const&) = delete;
        Generic& operator=(Generic&&) = delete;
        ~Generic() = default;
        
        bool apply(double const urn, mch::WangLandau<Value>& wangLandau, data::Data<Value> const& data, state::State<Value>& state, ut::UniformRng& urng, imp::itf::Batcher<Value>& batcher) {
            if(flag_ != ut::Flag::Pending) flag_ = prepare(urn, wangLandau, data, state, urng);
            if(flag_ == ut::Flag::Pending) flag_ = decide(wangLandau, data, state, batcher);
            return flag_ != ut::Flag::Pending;
        };
        
        jsx::value json() {
            return update_.json();
        };
        
    private:
        ut::Flag flag_;
        
        Upd update_;
        
        ImplementImpurity<Upd, Mode, Value> impurity_;
        ImplementBath<Upd, Value> bath_;
        
        ut::Zahl<double> x_;
        
        ut::Flag prepare(double const urn, mch::WangLandau<Value>& wangLandau, data::Data<Value> const& data, state::State<Value>& state, ut::UniformRng& urng) {
            if(update_.propose(urn, data, state, urng)) {
                if(impurity_.surviving(update_, data, state)) {
                    auto const ratio = (mch::itf::Update<Value>::ratioChoose() * update_.ratio(data, state) *
                                        impurity_.ratio(update_, data, state) * bath_.ratio(update_, data, state) *
                                        wangLandau.eta(update_.target(cfg::get<Origin>(state.worm()))) / wangLandau.eta(cfg::get<Origin>(state.worm())));
                    
                    if(!(ratio == .0)) {
                        x_ = urng()/ratio; return ut::Flag::Pending;
                    }
                    
                    bath_.reject(update_, data, state);
                }
                impurity_.reject(update_, data, state);
            }
            update_.reject(data, state);
            
            wangLandau.inc(cfg::get<Origin>(state.worm()));
            
            return ut::Flag::Reject;
        };
        
        ut::Flag decide(mch::WangLandau<Value>& wangLandau, data::Data<Value> const& data, state::State<Value>& state, imp::itf::Batcher<Value>& batcher) {
            auto const flag = impurity_.decide(x_, data, state, batcher);
            
            if(flag == ut::Flag::Reject) {
                update_.reject(data, state);
                impurity_.reject(update_, data, state);
                bath_.reject(update_, data, state);

                wangLandau.inc(cfg::get<Origin>(state.worm()));
            }
            
            if(flag == ut::Flag::Accept) {
                update_.accept(data, state);
                impurity_.accept(update_, data, state);
                bath_.accept(update_, data, state);
                
                wangLandau.inc(cfg::get<Target>(state.worm()));
            }
            
            return flag;
        };
    };
    
}

#endif
