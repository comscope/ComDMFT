#ifndef STATE_H
#define STATE_H

#include <vector>

#include "Utilities.h"
#include "impurity/Product.h"
#include "impurity/Dynamic.h"
#include "impurity/DensityMatrix.h"
#include "impurity/Exponentials.h"
#include "bath/Bath.h"
#include "config/Config.h"
#include "Data.h"


namespace state {
    
    struct State {
        State() = delete;
        template<typename Alloc, typename HybVal>
        State(jsx::value const& jParams, data::Data const& data, std::int64_t mcId, ut::Options<Alloc, HybVal>) :
        mcId_(mcId),
        config_(mcId_, data.ops().flavors()),
        product_(new imp::Product<Alloc>(jParams, data.eig(), data.ide(), data.ops())),
        densityMatrix_(new imp::DensityMatrix<Alloc>()),
        baths_(data.hyb().blocks().size()),
        dyn_(data.dyn() != nullptr ? new imp::Dynamic(*data.dyn()) : new imp::itf::Dynamic()),
        exponentials_(jParams, data){
            for(auto& bath : baths_)
                bath.reset(new bath::Bath<HybVal>());

            int index = 0;
            for(auto const& block : data.hyb().blocks()) {
                for(auto const& flavor : block.flavorsL())
                    for(auto const& entry : config_.at(flavor))
                        bath(index).insertL(entry.key(), flavor, entry.ptr());

                for(auto const& flavor : block.flavorsR())
                    for(auto const& entry : config_.at(flavor))
                        bath(index).insertR(entry.key(), flavor, entry.ptr());

                bath(index++).clean(data.hyb());
            }
            
            for(int b = 0; b < data.hyb().blocks().size(); ++b) {
                auto const& opsL = bath(b).opsL();
                auto const& opsR = bath(b).opsR();
                
                for(std::size_t i = 0; i < opsL.size(); ++i) {
                    if(!product().insert(opsR[i].key(), opsR[i].flavor(), opsR[i].ptr()))
                        throw std::runtime_error("state::State::constructor: key in config appears twice.");
                    if(!product().insert(opsL[i].key(), opsL[i].flavor(), opsL[i].ptr()))
                        throw std::runtime_error("state::State::constructor: key in config appears twice.");
                    
                    dyn().insert(opsR[i].key(), opsR[i].flavor()%2);
                    dyn().insert(opsL[i].key(), opsL[i].flavor()%2);
                }
            }
            
            dyn().ratio();
            dyn().accept();
        };
        State(State const&) = delete;
        State(State&&) = delete;
        State& operator=(State const&) = delete;
        State& operator=(State&&) = delete;
        ~State() = default;

    public:
        cfg::Config& config() { return config_;};
        cfg::Config const& config() const { return config_;};
        
        imp::itf::Product& product() { return *product_;};
        imp::itf::Product const& product() const { return *product_;};

        imp::itf::DensityMatrix& densityMatrix(){ return *densityMatrix_;};
        imp::itf::DensityMatrix const& densityMatrix() const { return *densityMatrix_;};

        bath::itf::Bath& bath(int i) { return *baths_[i];};
        bath::itf::Bath const& bath(int i) const { return *baths_[i];};
        
        imp::itf::Dynamic& dyn() { return *dyn_;};
        imp::itf::Dynamic const& dyn() const { return *dyn_;};
        
        void touch() {
            touched_ = true;
        };
        imp::Exponentials& exponentials() {
            if(touched_) { exponentials_.set(config_); touched_ = false; }
            return exponentials_;
        };

        ut::complex sign() const {
            ut::complex sign = product().perm()*(densityMatrix().Z() <= .0 ? -1 : 1);
            for(int b = 0; b < baths_.size(); ++b) sign *= bath(b).sign();
            return sign;
        };
        
        void clean(data::Data const& data) {
            dyn().clean();
            for(int b = 0; b < baths_.size(); ++b) bath(b).clean(data.hyb());
            exponentials_.clean(config_); 
        };

        std::int64_t mc_id() const {
            return mcId_;
        }
        
    private:
        std::int64_t const mcId_;

        cfg::Config config_;
        std::unique_ptr<imp::itf::Product> product_;
        std::unique_ptr<imp::itf::DensityMatrix> densityMatrix_;
        std::vector<std::unique_ptr<bath::itf::Bath>> baths_;
        std::unique_ptr<imp::itf::Dynamic> dyn_;
        imp::Exponentials exponentials_;
        
        bool touched_ = true;
    };
    
    
    namespace itf {
        
        struct Init {
            virtual bool apply(data::Data const&, State&, imp::itf::Batcher&) = 0;
            virtual ~Init() = default;
        };
        
    }
    
    template<typename Alloc>
    struct Init : itf::Init {
        Init() : flag_(ut::Flag::Reject) {};
        Init(Init const&) = delete;
        Init(Init&&) = delete;
        Init& operator=(Init const&) = delete;
        Init& operator=(Init&&) = delete;
        ~Init() = default;
        
        bool apply(data::Data const& data, State& state, imp::itf::Batcher& batcher) {
            if(flag_ != ut::Flag::Pending) flag_ = prepare(data, state);
            if(flag_ == ut::Flag::Pending) flag_ = decide(state, batcher);
            return flag_ != ut::Flag::Pending;
        };
        
    private:
        ut::Flag flag_;
        imp::DensityMatrix<Alloc> densityMatrix_;
        
        ut::Flag prepare(data::Data const& data, State& state) {
            densityMatrix_ = imp::DensityMatrix<Alloc>(state.product(), data.eig());
            if(ut::Flag::Pending != densityMatrix_.surviving(data.eig()))
                throw std::runtime_error("state::Init: initial trace is zero");
            return ut::Flag::Pending;
        };
        
        ut::Flag decide(State& state, imp::itf::Batcher& batcher) {
            if(ut::Flag::Pending == densityMatrix_.decide(.0, state.product(), batcher)) return ut::Flag::Pending;
            state.product().accept();
            imp::get<Alloc>(state.densityMatrix()) = std::move(densityMatrix_);
            return ut::Flag::Accept;
        };
    };
}

#endif
