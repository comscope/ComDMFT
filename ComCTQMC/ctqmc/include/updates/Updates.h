#ifndef UPDATES_IMPLEMENTATION_H
#define UPDATES_IMPLEMENTATION_H

#include "Config.h"
#include "Impurity.h"
#include "Bath.h"

#include "../Utilities.h"
#include "../Data.h"
#include "../State.h"

namespace upd {
    
    namespace itf {
        
        struct Update {
            virtual ~Update() = default;
        };
        
    };
    
    template<typename Def> struct Update;
    
    template<typename Def>
    struct Update : itf::Update {
        Update() = delete;
        template<typename Alloc, typename HybVal>
        Update(ut::Options<Alloc, HybVal>) : flag_(ut::Flag::Reject), confUpd_(new cfg::Update<Def>()), impUpd_(new imp::Update<Def, Alloc>()), bathUpd_(new bath::Update<Def, HybVal>()) {};
        Update(Update const&) = delete;
        Update(Update&&) = delete;
        Update& operator=(Update const&) = delete;
        Update& operator=(Update&&) = delete;
        ~Update() = default;
        
        bool operator()(Def& u, data::Data const& data, state::State& state, ut::UniformRng& urng, imp::itf::Batcher& batcher) {
            if(flag_ != ut::Flag::Pending) flag_ = prepare(u, data, state, urng);
            if(flag_ == ut::Flag::Pending) flag_ = decide(u, data, state, batcher);
            return flag_ != ut::Flag::Pending;
        };
        
    private:
        ut::Flag flag_;
        ut::Zahl x_;
        
        std::unique_ptr<cfg::Update<Def>> confUpd_;
        std::unique_ptr<imp::itf::Update<Def>> impUpd_;
        std::unique_ptr<bath::itf::Update<Def>> bathUpd_;
        
        ut::Flag prepare(Def& u, data::Data const& data, state::State& state, ut::UniformRng& urng) {
            if(confUpd_->propose(u, data, state, urng)) {
                if(ut::Flag::Pending == impUpd_->surviving(u, data, state)) {
                    auto const ratio = confUpd_->ratio(u, data, state) * bathUpd_->ratio(u, data, state) * impUpd_->ratio(u, data, state);
                    
                    if(!(ratio == .0)) {
                        x_ = urng()/ratio; return ut::Flag::Pending;
                    }

                    bathUpd_->reject(u, data, state);  //besser wenn alli updates da sind, au wenn's nuet machend
                }
                impUpd_->reject(u, data, state);
            }
            confUpd_->reject(u, data, state);
            return ut::Flag::Reject;
        };
        
        ut::Flag decide(Def const& u, data::Data const& data, state::State& state, imp::itf::Batcher& batcher) {
            auto const flag = impUpd_->decide(x_, data, state, batcher);
            
            if(flag == ut::Flag::Reject) {
                confUpd_->reject(u, data, state);
                impUpd_->reject(u, data, state);
                bathUpd_->reject(u, data, state);
            }
            
            if(flag == ut::Flag::Accept) {
                confUpd_->accept(u, data, state);
                impUpd_->accept(u, data, state);
                bathUpd_->accept(u, data, state);
            }
            
            return flag;
        };
    };
    
}

#endif
