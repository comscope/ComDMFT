#ifndef CTQMC_INCLUDE_MARKOVCHAIN_MARKOVCHAIN_H
#define CTQMC_INCLUDE_MARKOVCHAIN_MARKOVCHAIN_H

#include <vector>
#include <array>

#include "Update.h"
#include "../Utilities.h"
#include "../Data.h"
#include "../State.h"


namespace mch {
    
    template<typename Value>
    struct MarkovChain {
        MarkovChain() = delete;
        template<typename Mode>
        MarkovChain(jsx::value const& jParams, std::int64_t ID, Mode) :
        clean_(jParams.is("clean") ? jParams("clean").int64() : 10000), steps_(0),
        init_(new state::Init<Mode, Value>()),
        urng_(ut::Engine((jParams.is("seed") ? jParams("seed").int64() : 41085) + (jParams.is("seed increment") ? jParams("seed increment").int64() : 857)*ID), ut::UniformDistribution(.0, 1.)),
        update_(nullptr) {
        };
        MarkovChain(MarkovChain const&) = delete;
        MarkovChain(MarkovChain&&) = delete;
        MarkovChain& operator=(MarkovChain const&) = delete;
        MarkovChain& operator=(MarkovChain&&) = delete;
        ~MarkovChain() = default;
        
        
        void add(std::unique_ptr<itf::Update<Value>> update) {
            if(update->origin() != update->target())
                throw std::runtime_error("mc::MarkovChain::add");
            
            add_update(std::move(update));
        };
        
        void add(std::unique_ptr<itf::Update<Value>> updateAB, std::unique_ptr<itf::Update<Value>> updateBA) {
            if(updateAB->origin() != updateBA->target() || updateBA->origin() != updateAB->target())
                throw std::runtime_error("mc::MarkovChain::add");
            
            updateAB->ratioChoose_ = updateBA->probChoose()/updateAB->probChoose();  //not yet properly normalized
            updateBA->ratioChoose_ = updateAB->probChoose()/updateBA->probChoose();  //not yet properly normalized
            
            add_update(std::move(updateAB));
            add_update(std::move(updateBA));
        };
        
        void finalize(state::State<Value> const& state) {
            for(auto& updates : allUpdates_)
                for(auto& update : updates)
                    update->ratioChoose_ *= allDistrs_[update->origin()].back()/allDistrs_[update->target()].back();  // now they are properly normalized
            
            choose_update(state);
        };
        
        bool init(data::Data<Value> const& data, state::State<Value>& state, imp::itf::Batcher<Value>& batcher) {
            return init_->apply(data, state, batcher);
        };
        
        bool cycle(mch::WangLandau<Value>& wangLandau,  data::Data<Value> const& data, state::State<Value>& state, imp::itf::Batcher<Value>& batcher) {
            if(!update_->apply(urn_, wangLandau, data, state, urng_, batcher)) return false;
            
            choose_update(state); if(++steps_ % clean_ == 0) state.clean(data);
            
            return true;
        };
        
    private:
        std::int64_t const clean_;
        std::int64_t steps_;
        double urn_;
        
        std::unique_ptr<state::itf::Init<Value>> init_;
        
        ut::UniformRng urng_;
        std::array<std::vector<double>, cfg::Worm::size()> allDistrs_;
        std::array<std::vector<std::unique_ptr<itf::Update<Value>>>, cfg::Worm::size()> allUpdates_;
        itf::Update<Value>* update_;

        
        void add_update(std::unique_ptr<itf::Update<Value>> update) {
            auto& distr = allDistrs_[update->origin()];
            auto& updates = allUpdates_[update->origin()];
            
            distr.push_back((distr.size() ? distr.back() : .0) + update->probChoose());
            updates.push_back(std::move(update));
        };
        
        void choose_update(state::State<Value> const& state) {
            auto& distr = allDistrs_[state.worm().index()];
            auto& updates = allUpdates_[state.worm().index()];
            
            double const urn = urng_()*distr.back();
            auto const it = std::upper_bound(distr.begin(), distr.end(), urn);
            double const low = it != distr.begin() ? *(it - 1) : .0;
            
            urn_ = (urn - low)/(*it - low);
            update_ = updates[it - distr.begin()].get();
        };
    };
    
}

#endif
