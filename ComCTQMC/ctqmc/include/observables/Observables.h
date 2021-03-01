#ifndef CTQMC_INCLUDE_OBSERVABLES_OBSERVABLES_H
#define CTQMC_INCLUDE_OBSERVABLES_OBSERVABLES_H

#include <vector>

#include "Observable.h"
#include "../Data.h"
#include "../State.h"

namespace obs {

    template<typename Value>
    struct WormObservables {
        WormObservables() = delete;
        WormObservables(std::string worm) :
        worm_(worm), steps_(0), it_(obs_.end()) {
        }
        WormObservables(WormObservables const&) = delete;
        WormObservables(WormObservables&&) = default;
        WormObservables& operator=(WormObservables const&) = delete;
        WormObservables& operator=(WormObservables&&) = delete;
        ~WormObservables() = default;
        
        
        template<typename T, typename... Args>
        void add(std::int64_t sweep, std::int64_t store, Args&&... args) {
            obs_.emplace_back(sweep, obs_pointer(new T(store, std::forward<Args>(args)...)));
            
            it_ = obs_.end();
        };
        
        bool sample(data::Data<Value> const& data, state::State<Value>& state) {
            if(it_ != obs_.end()) return false;
            
            ++steps_;
            
            for(auto const& obs : obs_)
                if(steps_%obs.first == 0)
                {
                    sign_ = state.sign();
                    it_   = obs_.begin();
                    
                    return true;
                }
            
            return false;
        }
        
        bool cycle(data::Data<Value> const& data, state::State<Value>& state, jsx::value& measurements, imp::itf::Batcher<Value>& batcher) {
            while(it_ != obs_.end())
            {
                if(steps_%it_->first == 0)
                    if(!it_->second->sample(sign_, data, state, measurements[worm_], batcher))
                        return false;

                ++it_;
            }
            
            return true;
        };
        
        void finalize(data::Data<Value> const& data, jsx::value& measurements) {
            for(auto& obs : obs_)
                obs.second->finalize(data, measurements[worm_]);
        };
        
    private:
        using obs_pointer = std::unique_ptr<itf::Observable<Value>>;
        
        std::string  const worm_;
        std::int64_t steps_;
        
        Value sign_;
        std::vector<std::pair<std::int64_t, obs_pointer>> obs_;
        typename std::vector<std::pair<std::int64_t, obs_pointer>>::iterator it_;
    };
    
    
    template<typename Value> using Observables = std::array<std::unique_ptr<WormObservables<Value>>, cfg::Worm::size()>;

}

#endif
