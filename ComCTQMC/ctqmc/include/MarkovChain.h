#ifndef MARKOVCHAIN
#define MARKOVCHAIN

#include <vector>

#include "Utilities.h"
#include "Updates.h"
#include "Data.h"
#include "Weight.h"
#include "Observables.h"
#include "Config.h"

// Irgendwas stimmt da no nit, MarkovChain is nit dr richtig name: MarkovChain -> MonteCarlo,  MonteCarlo -> Simulation
// isch doch alles hippie-kacke hippie-kacke hippie-kacke

namespace ma {
    
    template<class GreenMeas, class HybFunc>
    struct MarkovChain {
        MarkovChain() = delete;
        MarkovChain(jsx::value const& jParams, int const id) :
        id_(id),
        urng_(ut::Engine((jParams.is("seed") ? jParams("seed").int64() : 41085) + (jParams.is("seed increment") ? jParams("seed increment").int64() : 857)*id_), ut::UniformDistribution(.0, 1.)),
        data_(*da::Data<HybFunc>::Instance(jParams)),
        weight_(*we::Weight<HybFunc>::Instance(jParams, data_, id)), //Soetti z'Config file uebercho ....
        obs_(*ob::Observables<GreenMeas, HybFunc>::Instance(jParams, data_, weight_)),
        config_(*new co::ConfigCSQ(jParams("hybridisation")("matrix").size(), weight_.baths())) {
            auto const& jMatrix = jParams("hybridisation")("matrix");
            
            for(std::size_t i = 0; i < jMatrix.size(); ++i) 
                for(std::size_t j = 0; j < jMatrix.size(); ++j)
                    if(jMatrix(i)(j).string() != "") {
                        update_.emplace_back(new Visitor<up::InsertTwo>(2*i + 1, 2*j, data_.hyb().block(i)));
                        prob_.push_back((prob_.size() ? prob_.back() : .0) + 1.);
                        
                        update_.emplace_back(new Visitor<up::EraseTwo>(2*i + 1, 2*j, data_.hyb().block(i)));
                        prob_.push_back((prob_.size() ? prob_.back() : .0) + 1.);
                    }
            
            update_.emplace_back(new InitVisitor());
            upd_ = update_.back().get();
        };
        MarkovChain(MarkovChain const&) = delete;
        MarkovChain(MarkovChain&&) = delete;
        MarkovChain& operator=(MarkovChain const&) = delete;
        MarkovChain& operator=(MarkovChain&&) = delete;
        
        int id() {
            return id_;
        }
        
        int udpate() {
            if(upd_->apply(*this) == tr::Flag::Pending) return 0;
            
            upd_ = update_[std::upper_bound(prob_.begin(), prob_.end(), urng_()*prob_.back()) - prob_.begin()].get();
            
            return 1;
        }
        
        void sample() {
            obs_.sample(data_, weight_);
        };
        void store(jsx::value& measurements) {
            obs_.store(measurements);
        };
        
        int csample() {
            return obs_.csample(data_, weight_);
        };        
        void cstore(jsx::value& measurements) {
            obs_.cstore(measurements, data_);
        };
        
        void clean() {
            weight_.clean(data_);
            obs_.clean();
        };
        
        ~MarkovChain() {
            delete &config_;
            
            ob::Observables<GreenMeas, HybFunc>::Destroy(&obs_);
            we::Weight<HybFunc>::Destroy(&weight_);
            da::Data<HybFunc>::Destroy(&data_);
        };
        
    private:
        
        struct AbstractVisitor {
            virtual tr::Flag apply(MarkovChain&) = 0; // Es wuerd gnueege wenn de kack da 0 oder 1 zrugg git ......
            virtual ~AbstractVisitor() {};
        };
        
        template<typename U>
        struct Visitor : AbstractVisitor {
            template<typename... Args>
            Visitor(Args... args) : flag_(tr::Flag::Reject), u_(args...) {};
            tr::Flag apply(MarkovChain& maCh) {
                if(flag_ != tr::Flag::Pending) flag_ = maCh.prepare(u_);
                if(flag_ == tr::Flag::Pending) flag_ = maCh.decide(u_);
                return flag_;
            }
        private:
            tr::Flag flag_;
            U u_;
        };
        
        struct InitVisitor : AbstractVisitor {    // Isch das die richtig loesig fuer de init scheiss ?
            InitVisitor() : flag_(tr::Flag::Reject) {};
            tr::Flag apply(MarkovChain& maCh) {
                if(flag_ != tr::Flag::Pending) flag_ = maCh.weight_.prepare();
                if(flag_ == tr::Flag::Pending) flag_ = maCh.weight_.decide();
                return flag_;
            }
        private:
            tr::Flag flag_;
        };
        
        int const id_;
        
        ut::UniformRng urng_;
        
        da::Data<HybFunc> const& data_;                 //delete
        we::Weight<HybFunc>& weight_;                   //delete
        ob::Observables<GreenMeas, HybFunc>& obs_;      //delete
        
        co::Abstract& config_;                          //delete
        
        up::Tuple<tr::Updates, up::InsertTwo, up::EraseTwo> traceUpds_;
        up::Tuple<ba::Updates, up::InsertTwo, up::EraseTwo> bathUpds_;
        up::Tuple<dy::Updates, up::InsertTwo, up::EraseTwo> dynUpds_;
        
        std::vector<double> prob_;
        std::vector<std::unique_ptr<AbstractVisitor>> update_;
        AbstractVisitor* upd_;
        
        // One could also make a class for each update, but not necessary (for now at least)
        template<typename U>
        tr::Flag prepare(U& u) {
            auto& traceUpd = up::get<U>(traceUpds_);
            auto& bathUpd = up::get<U>(bathUpds_);
            auto& dynUpd = up::get<U>(dynUpds_);
            
            if(config_.propose(u, urng_)) {
                if(tr::Flag::Pending == traceUpd.surviving(u, weight_.trace())) {
                    za::Zahl const ratio = config_.ratio(u) * bathUpd.ratio(u, weight_.baths(), data_.hyb()) * dynUpd.ratio(u, weight_.dyn());

                    if(!(ratio == .0)) {
                        traceUpd.prepare(weight_.trace(), urng_()/ratio);
                        return tr::Flag::Pending;
                    }

                    bathUpd.reject(u, weight_.baths());  //besser wenn alli updates da sind, au wenn's nuet machend
                    dynUpd.reject(weight_.dyn());
                }
                traceUpd.reject(weight_.trace());
            }
            
            config_.reject(u);
            return tr::Flag::Reject;
        };
        
        template<typename U>
        tr::Flag decide(U const& u) {
            auto& traceUpd = up::get<U>(traceUpds_);
            auto& bathUpd = up::get<U>(bathUpds_);
            auto& dynUpd = up::get<U>(dynUpds_);
            
            tr::Flag const flag = traceUpd.decide(weight_.trace());
            
            if(flag == tr::Flag::Reject) {
                config_.reject(u);
                traceUpd.reject(weight_.trace());
                bathUpd.reject(u, weight_.baths());
                dynUpd.reject(weight_.dyn());
            }
            
            if(flag == tr::Flag::Accept) {
                config_.accept(u);
                traceUpd.accept(weight_.trace());
                bathUpd.accept(u, weight_.baths());
                dynUpd.accept(weight_.dyn());
                
                obs_.update(u);
            }
            
            return flag;
        };
    };
}

#endif
