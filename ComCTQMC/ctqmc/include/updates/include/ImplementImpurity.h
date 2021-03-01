#ifndef CTQMC_INCLUDE_UPDATES_INCLUDE_IMPLEMENTIMPURITY_H
#define CTQMC_INCLUDE_UPDATES_INCLUDE_IMPLEMENTIMPURITY_H


#include "../../Data.h"
#include "../../State.h"


namespace upd {
    
    template<typename Update, typename Mode, typename Value, typename = void>
    struct ImplementImpurity {
        bool surviving(Update const& update, data::Data<Value> const& data, state::State<Value>& state) const {
            return true;
        };
        
        ut::Zahl<double> ratio(Update const& update, data::Data<Value> const& data, state::State<Value>& state) const {
            return 1.;
        };
        
        ut::Flag decide(ut::Zahl<double> const& x, data::Data<Value> const& data, state::State<Value>& state, imp::itf::Batcher<Value>& batcher) const {
            return 1. <= x ? ut::Flag::Reject : ut::Flag::Accept;
        };
        
        void accept(Update const& update, data::Data<Value> const& data, state::State<Value>& state) const {
        };
        
        void reject(Update const& update, data::Data<Value> const& data, state::State<Value>& state) const {
        };
    };
    
    
    template<typename Update, typename Mode, typename Value>
    struct ImplementImpurity<Update, Mode, Value, ut::void_t<decltype(&Update::template impurity<Mode, Value>)>> {
        bool surviving(Update const& update, data::Data<Value> const& data, state::State<Value>& state) {
            if(!update.template impurity<Mode>(data, state)) return false;
            
            densityMatrix_ = imp::DensityMatrix<Mode, Value>(state.product(), data.eig());
            return ut::Flag::Pending == densityMatrix_.surviving(data.eig());
        };
        
        ut::Zahl<double> ratio(Update const& update, data::Data<Value> const& data, state::State<Value>& state) const {
            return state.dyn().ratio()*state.fact().ratio();
        };
        
        ut::Flag decide(ut::Zahl<double> const& x, data::Data<Value> const& data, state::State<Value>& state, imp::itf::Batcher<Value>& batcher) {
            return densityMatrix_.decide(ut::abs(state.densityMatrix().Z())*x, state.product(), batcher);
        };
        
        void accept(Update const& update, data::Data<Value> const& data, state::State<Value>& state) {
            imp::get<Mode>(state.densityMatrix()) = std::move(densityMatrix_);
            state.signTimeOrder() *= state.product().accept();
            state.fact().accept();
            state.dyn().accept();
        };
        
        void reject(Update const& update, data::Data<Value> const& data, state::State<Value>& state) {
            densityMatrix_ = imp::DensityMatrix<Mode, Value>();
            state.product().reject();
            state.fact().reject();
            state.dyn().reject();
        };
    private:
        imp::DensityMatrix<Mode, Value> densityMatrix_;
    };
    
}

#endif
