#ifndef CTQMC_INCLUDE_STATE_H
#define CTQMC_INCLUDE_STATE_H

#include <vector>

#include "Utilities.h"
#include "Data.h"
#include "impurity/Product.h"
#include "impurity/Dynamic.h"
#include "impurity/DensityMatrix.h"
#include "impurity/Fact.h"
#include "bath/Bath.h"
#include "config/Expansion.h"
#include "config/Worms.h"


namespace state {
    
    template<typename Mode>
    struct insert_worm_functor {
        
        template<typename Worm, typename Value>
        void operator()(Worm const& worm, data::Data<Value> const& data, state::State<Value>& state) {
            if(!cfg::insert_worm<Mode>(worm, data, state))
                throw std::runtime_error("state::insert_worm_functor: key appears twice in " + Worm::name() + " worm");
        };

    };
    
    
    template<typename Value>
    struct State {
        State() = delete;
        template<typename Mode>
        State(jsx::value const& jParams, data::Data<Value> const& data, jsx::value& jConfig, Mode) :
        signTimeOrder_(1),
        expansion_(jConfig.is("expansion") ? jConfig("expansion") : jsx::null_t(), data.ops().flavors()),
        worm_(jConfig.is("worm") ? jConfig("worm") : jsx::object_t{{"name", "partition"}, {"entry", jsx::null_t()}}),
        product_(new imp::Product<Mode, Value>(jParams, data.eig(), data.ide(), data.ops())),
        densityMatrix_(new imp::DensityMatrix<Mode, Value>()),
        baths_(data.hyb().blocks().size()),
        dyn_(data.dyn() != nullptr ? new imp::Dynamic(*data.dyn()) : new imp::itf::Dynamic()) {
            int index = 0;
            for(auto const& block : data.hyb().blocks()) {
                for(auto const& flavor : block.flavorsL())
                    for(auto const& key : expansion_.at(flavor))
                        baths_[index].insertL(key, flavor);
                
                for(auto const& flavor : block.flavorsR())
                    for(auto const& key : expansion_.at(flavor))
                        baths_[index].insertR(key, flavor);
                
                baths_[index++].clean(data.hyb());
            }
            
            for(auto const& bath : baths()) {
                auto const& opsL = bath.opsL();
                auto const& opsR = bath.opsR();
                
                for(std::size_t i = 0; i < opsL.size(); ++i) {
                    if(!product().insert(opsR[i].key(), opsR[i].flavor()))
                        throw std::runtime_error("state::State::constructor: key in config appears twice.");
                    if(!product().insert(opsL[i].key(), opsL[i].flavor()))
                        throw std::runtime_error("state::State::constructor: key in config appears twice.");
                    
                    dyn().insert(opsR[i].key(), opsR[i].flavor());
                    dyn().insert(opsL[i].key(), opsL[i].flavor());
                }
            }

            cfg::apply(worm(), insert_worm_functor<Mode>(), data, *this);


            dyn().ratio();
            dyn().accept();
            
            fact().accept();
        };
        State(State const&) = delete;
        State(State&&) = delete;
        State& operator=(State const&) = delete;
        State& operator=(State&&) = delete;
        ~State() = default;
        
    public:
        
        int& signTimeOrder() {
            return signTimeOrder_;
        };
        cfg::Expansion& expansion() {
            return expansion_;
        };
        cfg::Worm& worm() {
            return worm_;
        };
        imp::itf::Product<Value>& product() {
            return *product_;
        };
        imp::itf::DensityMatrix<Value>& densityMatrix() {
            return *densityMatrix_;
        };
        std::vector<bath::Bath<Value>>& baths() {
            return baths_;
        };
        imp::itf::Dynamic& dyn() {
            return *dyn_;
        };
        imp::Fact<Value>& fact() {
            return fact_;
        };
        
        
        int signTimeOrder() const {
            return signTimeOrder_;
        };
        cfg::Expansion const& expansion() const {
            return expansion_;
        };
        cfg::Worm const& worm() const {
            return worm_;
        };
        imp::itf::Product<Value> const& product() const {
            return *product_;
        };
        imp::itf::DensityMatrix<Value> const& densityMatrix() const {
            return *densityMatrix_;
        };
        std::vector<bath::Bath<Value>> const& baths() const {
            return baths_;
        };
        imp::itf::Dynamic const& dyn() const {
            return *dyn_;
        };
        imp::Fact<Value> const& fact() const {
            return fact_;
        };
        
        
        Value sign() const {
            Value sign = static_cast<double>(signTimeOrder())*densityMatrix().sign()*fact().sign();  //!!!!!!!!!!!!!!!!
            for(auto const& bath : baths()) sign *= bath.sign();
            return sign;
        };
        
        void clean(data::Data<Value> const& data) {
            dyn().clean();
            for(auto& bath : baths())
                bath.clean(data.hyb());
        };
        
        jsx::value json() const {
            return jsx::object_t{
                {"expansion", expansion().json()},
                {"worm", worm().json()}
            };
        };
        
    private:
        int signTimeOrder_;
        cfg::Expansion expansion_;
        cfg::Worm worm_;
        
        std::unique_ptr<imp::itf::Product<Value>> product_;
        std::unique_ptr<imp::itf::DensityMatrix<Value>> densityMatrix_;
        std::vector<bath::Bath<Value>> baths_;
        std::unique_ptr<imp::itf::Dynamic> dyn_;
        imp::Fact<Value> fact_;
    };
    
    
    
    namespace itf {
        
        template<typename Value>
        struct Init {
            virtual bool apply(data::Data<Value> const&, State<Value>&, imp::itf::Batcher<Value>&) = 0;
            virtual ~Init() = default;
        };
        
    }
    
    template<typename Mode, typename Value>
    struct Init : itf::Init<Value> {
        Init() : flag_(ut::Flag::Reject) {};
        Init(Init const&) = delete;
        Init(Init&&) = delete;
        Init& operator=(Init const&) = delete;
        Init& operator=(Init&&) = delete;
        ~Init() = default;
        
        bool apply(data::Data<Value> const& data, State<Value>& state, imp::itf::Batcher<Value>& batcher) {
            if(flag_ != ut::Flag::Pending) flag_ = prepare(data, state);
            if(flag_ == ut::Flag::Pending) flag_ = decide(state, batcher);
            return flag_ != ut::Flag::Pending;
        };
        
    private:
        ut::Flag flag_;
        imp::DensityMatrix<Mode, Value> densityMatrix_;
        
        ut::Flag prepare(data::Data<Value> const& data, State<Value>& state) {
            densityMatrix_ = imp::DensityMatrix<Mode, Value>(state.product(), data.eig());
            if(ut::Flag::Pending != densityMatrix_.surviving(data.eig()))
                throw std::runtime_error("state::Init: initial trace is zero");
            return ut::Flag::Pending;
        };
        
        ut::Flag decide(State<Value>& state, imp::itf::Batcher<Value>& batcher) {
            if(ut::Flag::Pending == densityMatrix_.decide(.0, state.product(), batcher)) return ut::Flag::Pending;
            imp::get<Mode>(state.densityMatrix()) = std::move(densityMatrix_);
            state.signTimeOrder() *= state.product().accept();
            
            return ut::Flag::Accept;
        };
    };
}

#endif
