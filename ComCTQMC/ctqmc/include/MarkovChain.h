#ifndef MARKOVCHAIN_H
#define MARKOVCHAIN_H

#include <vector>

#include "Utilities.h"
#include "Data.h"
#include "State.h"
#include "updates/Definition.h"
#include "updates/Updates.h"


namespace mc {
    
    struct TypeId {
        using FuncPtr = TypeId(*)();
        TypeId(FuncPtr val) : val_(val) {};
    private:
        FuncPtr val_;
        
        friend bool operator<(TypeId const& lhs, TypeId const& rhs) {
            return lhs.val_ < rhs.val_;
        };
    };
    
    template<typename T> TypeId type_id() {
        return &type_id<T>;
    };
    
    
    struct MarkovChain {
        MarkovChain() = delete;
        MarkovChain(jsx::value const& jParams, std::int64_t mcId) :
        clean_(jParams.is("clean") ? jParams("clean").int64() : 10000),
        steps_(0),
        urng_(ut::Engine((jParams.is("seed") ? jParams("seed").int64() : 41085) + (jParams.is("seed increment") ? jParams("seed increment").int64() : 857)*mcId), ut::UniformDistribution(.0, 1.)),
        upd_(nullptr) {
        };
        MarkovChain(MarkovChain const&) = delete;
        MarkovChain(MarkovChain&&) = delete;
        MarkovChain& operator=(MarkovChain const&) = delete;
        MarkovChain& operator=(MarkovChain&&) = delete;
        ~MarkovChain() = default;
        
        
        template<typename Alloc, typename HybVal, typename Def>
        void add(Def def, double prob) {
            updates_.push_back(make_update<Alloc, HybVal>(def));
            
            prob_.push_back((prob_.size() ? prob_.back() : .0) + prob);
        };
        
        
        bool cycle(data::Data const& data, state::State& state, imp::itf::Batcher& batcher) {
            if(upd_ == nullptr) upd_ = updates_[std::upper_bound(prob_.begin(), prob_.end(), urng_()*prob_.back()) - prob_.begin()].get();
            
            if(!upd_->apply(data, state, urng_, batcher)) return false;
            
            upd_ = nullptr; ++steps_; if(steps_ % clean_ == 0) state.clean(data);
            
            return true;
        };
        
        
        std::int64_t steps() const {
            return steps_;
        };
        
    private:
        
        struct BindItf {
            virtual bool apply(data::Data const&, state::State&, ut::UniformRng&, imp::itf::Batcher&) = 0;
            virtual ~BindItf() = default;
        };
        
        template<typename Def>
        struct Bind : BindItf {
            Bind(Def def, upd::Update<Def>& impl) : def_(def), impl_(impl) {};
            
            bool apply(data::Data const& data, state::State& state, ut::UniformRng& urng, imp::itf::Batcher& batcher) {
                return impl_(def_, data, state, urng, batcher);
            };
        private:
            Def def_;
            upd::Update<Def>& impl_;
        };
        
        std::int64_t const clean_;
        std::int64_t steps_;
        
        ut::UniformRng urng_;
        std::vector<double> prob_;
        std::vector<std::unique_ptr<BindItf>> updates_;
        BindItf* upd_;
        
        std::map<TypeId, std::unique_ptr<upd::itf::Update>> impl_;
       
        
        template<typename Alloc, typename HybVal, typename Def>
        std::unique_ptr<BindItf> make_update(Def def) {
            using Impl = upd::Update<Def>; auto const typeId = type_id<Impl>();
            
            if(!impl_.count(typeId)) impl_[typeId].reset(new Impl(ut::Options<Alloc, HybVal>()));

            return std::unique_ptr<BindItf>(new Bind<Def>(def, static_cast<Impl&>(*impl_[typeId])));
        }
    };
    
    
    
    
    template<typename Alloc, typename HybVal>
    void addUpdates(jsx::value const& jParams, data::Data const& data, MarkovChain& markovChain) {
        mpi::cout << "Setting default updates ... ";
        
        int bath = 0;
        for(auto const& block : data.hyb().blocks()) {
            for(auto const& flavorL : block.flavorsL())
                for(auto const& flavorR : block.flavorsR()) {
                    markovChain.add<Alloc, HybVal>(upd::InsertTwo(flavorL, flavorR, bath), 1.);
                    markovChain.add<Alloc, HybVal>(upd::EraseTwo(flavorL, flavorR, bath), 1.);
                }
            ++bath;
        }
        
        mpi::cout << "Ok" << std::endl;
    };
    
}

#endif
