#ifndef CTQMC_INCLUDE_MARKOVCHAIN_UPDATE_H
#define CTQMC_INCLUDE_MARKOVCHAIN_UPDATE_H

#include "WangLandau.h" 
#include "../Utilities.h"
#include "../Data.h"
#include "../State.h"
#include "../../../include/JsonX.h"

namespace mch {
    
    template<typename> struct MarkovChain;
    
    namespace itf {
        
        template<typename Value>
        struct Update {
            Update() = delete;
            Update(int origin, int target, double probChoose) : origin_(origin), target_(target), probChoose_(probChoose), ratioChoose_(1.) {};
            Update(Update const&) = delete;
            Update(Update&&) = delete;
            Update& operator=(Update const&) = delete;
            Update& operator=(Update&&) = delete;
            virtual ~Update() = default;
            
            int origin() const { return origin_;};
            int target() const { return target_;};
            
            double probChoose() const { return probChoose_;};
            double ratioChoose() const { return ratioChoose_;};
            
            virtual bool apply(double const, mch::WangLandau<Value>&, data::Data<Value> const&, state::State<Value>&, ut::UniformRng&, imp::itf::Batcher<Value>&) = 0;

            virtual jsx::value json() = 0;  // should return statistics about move acceptance/rejectance
        
        private:
            int const origin_, target_;
            double const probChoose_;
            double ratioChoose_;
            
            template<typename> friend struct mch::MarkovChain;
        };
        
    };
    
    template<typename Value> using unique_update_ptr = std::unique_ptr<mch::itf::Update<Value>>;
    
}

#endif
