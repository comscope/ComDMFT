#ifndef CTQMC_INCLUDE_OBSERVABLES_PARTITION_SETUP_H
#define CTQMC_INCLUDE_OBSERVABLES_PARTITION_SETUP_H

#include <vector>

#include "OneParticle.h"
#include "Misc.h"
#include "SectorProb.h"
#include "SuscFlavor.h"
#include "Chain.h"
#include "DensityMatrix.h"
#include "SuscMatrix.h"

#include "../Observables.h"

namespace obs {
    
    namespace partition {
        
        template<typename Mode, typename Value>
        void setup_obsA(jsx::value const& jParams, data::Data<Value>& data, WormObservables<Value>& observables) {
            jsx::value jPartition = jParams(cfg::partition::Worm::name());
            
            
            std::int64_t const sweepA = jPartition.is("sweepA") ? jPartition("sweepA").int64() : 50;
            std::int64_t const storeA = jPartition.is("storeA") ? jPartition("storeA").int64() : 100;

            
            observables.template add<Misc<Value>>(sweepA, storeA, jParams, data);
            
            if(jPartition.is("green basis") ? jPartition("green basis").string() == "matsubara" : true)
                observables.template add<OneParticle<BinMoments, Green, Value>>(sweepA, storeA, jParams, data);
            else if(jPartition("green basis").string() == "legendre")
                observables.template add<OneParticle<Legendre, Green, Value>>(sweepA, storeA, jParams, data);
            else
                throw std::runtime_error("addObservablesA: unknown green basis option");
            
            if( ! (jPartition.is("density matrix precise") ? jPartition("density matrix precise").boolean() : false))
                observables.template add<DensityMatrixStatic<Mode, Value>>(sweepA, storeA, jParams, data);
            
            if(data.dyn() != nullptr) observables.template add<DensityMatrixDynamic<Mode, Value>>(sweepA, storeA, jParams, data);
            

            if((jPartition.is("occupation susceptibility direct") ? jPartition("occupation susceptibility direct").boolean() : false) && 
               (jPartition.is("occupation susceptibility bulla")  ? jPartition("occupation susceptibility bulla").boolean()  : false))
                observables.template add<SuscMatrix<Mode, Value, true, true>>(sweepA, storeA, jParams, data);
            else  if(jPartition.is("occupation susceptibility direct") ? jPartition("occupation susceptibility direct").boolean() : false)
                observables.template add<SuscMatrix<Mode, Value, true, false>>(sweepA, storeA, jParams, data);
            else if(jPartition.is("occupation susceptibility bulla")  ? jPartition("occupation susceptibility bulla").boolean()  : false)
                observables.template add<SuscMatrix<Mode, Value, false, true>>(sweepA, storeA, jParams, data);
        }
        
        
        template<typename Mode, typename Value>
        void setup_obsB(jsx::value const& jParams, data::Data<Value>& data, WormObservables<Value>& observables) {
            jsx::value jPartition = jParams(cfg::partition::Worm::name());

            std::int64_t const sweepB = jPartition.is("sweepB") ? jPartition("sweepB").int64() : 250;
            std::int64_t const storeB = jPartition.is("storeB") ? jPartition("storeB").int64() : 20;
            
            observables.template add<SectorProb<Mode, Value>>(sweepB, storeB, jParams, data);
            
            if((jPartition.is("quantum number susceptibility") ? jPartition("quantum number susceptibility").boolean() : false) ||
               (jPartition.is("occupation susceptibility bulla")  ? jPartition("occupation susceptibility bulla").boolean()  : false))
                observables.template add<SuscFlavor<Value>>(sweepB, storeB, jParams, data);
            
            if((jPartition.is("density matrix precise") ? jPartition("density matrix precise").boolean() : false) &&
               (jPartition.is("green bulla") ? jPartition("green bulla").boolean() : true))
                observables.template add<Chain<Mode, Value, true, true>>(sweepB, storeB, jParams, data);
            else if(jPartition.is("density matrix precise") ? jPartition("density matrix precise").boolean() : false)
                observables.template add<Chain<Mode, Value, true, false>>(sweepB, storeB, jParams, data);
            else if(jPartition.is("green bulla") ? jPartition("green bulla").boolean() : true)
                observables.template add<Chain<Mode, Value, false, true>>(sweepB, storeB, jParams, data);
            
            if(jPartition.is("green bulla") ? jPartition("green bulla").boolean() : true) {
                if(jPartition.is("green basis") ? jPartition("green basis").string() == "matsubara" : true) {
                    observables.template add<OneParticle<BinMoments, BullaL, Value>>(sweepB, storeB, jParams, data);
                    observables.template add<OneParticle<BinMoments, BullaR, Value>>(sweepB, storeB, jParams, data);
                } else if(jPartition("green basis").string() == "legendre") {
                    observables.template add<OneParticle<Legendre, BullaL, Value>>(sweepB, storeB, jParams, data);
                    observables.template add<OneParticle<Legendre, BullaR, Value>>(sweepB, storeB, jParams, data);
                } else
                    throw std::runtime_error("obs::Observables: unknown green basis option");
            }
        }
        
    }
    
    
    template<typename Mode, typename Value>
    void setup_worm_obs(jsx::value const& jParams, data::Data<Value>& data, std::unique_ptr<WormObservables<Value>>& observables, ut::wrap<cfg::partition::Worm>) {
        jsx::value jPartition = jParams(cfg::partition::Worm::name());

        
        if(jPartition.is("green bulla") ? jPartition("green bulla").boolean() : true) {
            auto& bullaOps = data.template opt<imp::itf::BullaOperators<Value>>();
            if(bullaOps.get() == nullptr) bullaOps.reset(new imp::BullaOperators<Mode, Value>(jParams("hloc")("interaction"), jParams("operators"), data.eig()));
        }
        
        if(jPartition.is("occupation susceptibility direct") ? jPartition("occupation susceptibility direct").boolean() : false) {
            auto& occ = data.template opt<imp::itf::Occupation<Value>>();
            if(occ.get() == nullptr) occ.reset(new imp::Occupation<Mode, Value>(jParams("operators"), data.eig()));
        }
        
        if(jPartition.is("occupation susceptibility bulla")  ? jPartition("occupation susceptibility bulla").boolean()  : false) {
            auto& bullaOcc = data.template opt<imp::itf::BullaOccupation<Value>>();
            bullaOcc.reset(new imp::BullaOccupation<Mode, Value>(jParams, data.filling(), jParams("hloc")("eigen values"), jParams("operators"), data.eig()));
        }
        
        
        observables.reset(new WormObservables<Value>(cfg::partition::Worm::name()));

        
        partition::setup_obsA<Mode>(jParams, data, *observables);
        partition::setup_obsB<Mode>(jParams, data, *observables);
    }
}

#endif
