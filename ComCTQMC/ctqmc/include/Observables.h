#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <vector>

#include "observables/OneParticle.h"
#include "observables/Misc.h"
#include "observables/SectorProb.h"
#include "observables/SuscMatrix.h"
#include "observables/SuscFlavor.h"
#include "observables/Chain.h"
#include "observables/DensityMatrix.h"

namespace obs {

    struct Observables {
        Observables() = delete;
        Observables(jsx::value const& jParams, std::int64_t store) :
        store_(store), samples_(0), lock_(false), phase_(Phase::Prepare) {
        }
        Observables(Observables const&) = delete;
        Observables(Observables&&) = delete;
        Observables& operator=(Observables const&) = delete;
        Observables& operator=(Observables&&) = delete;
        ~Observables() = default;
        
        bool& lock() { return lock_;};
        
        template<typename T, typename... Args>
        void add(Args&&... args) {
            obs_.emplace_back(new T(std::forward<Args>(args)...));
        };
        
        bool sample(data::Data const& data, state::State& state, jsx::value& measurements, imp::itf::Batcher& batcher) {
            switch(phase_) {
                case Phase::Prepare:
                    state.touch();
                    
                    sign_ = state.sign(); it_ = obs_.begin();
                    
                    phase_ = Phase::Sample;
                    
                case Phase::Sample:
                    while(it_ != obs_.end()) {
                        if(!(*it_)->sample(sign_, data, state, batcher)) return false;
                        ++it_;
                    }
                    ++samples_;
                    
                    if(samples_%store_ == 0) {
                        for(auto& obs : obs_) obs->store(data, measurements, samples_);
                        samples_ = 0;
                    }
                    phase_ = Phase::Prepare;
            }
            
            return true;
        };
        
        void finalize(data::Data const& data, jsx::value& measurements) {
            for(auto& obs : obs_) obs->finalize(data, measurements, samples_);
            samples_ = 0;
        };
        
    private:
        
        enum class Phase { Prepare, Sample };
        
        std::int64_t const store_;
        std::int64_t samples_;
        
        bool lock_;
        Phase phase_;
        
        ut::complex sign_;
        std::vector<std::unique_ptr<Observable>> obs_;
        std::vector<std::unique_ptr<Observable>>::iterator it_;
    };
    
    
    
    template<typename Alloc, typename HybVal>
    void addObservablesA(jsx::value const& jParams, data::Data const& data, Observables& obs) {
        using namespace obs;
        
        obs.add<Misc>(jParams, data);
        
        if(jParams.is("green basis") ? jParams("green basis").string() == "matsubara" : true)
            obs.add<OneParticle<BinMoments, Green, HybVal>>(jParams, data);
        else if(jParams("green basis").string() == "legendre")
            obs.add<OneParticle<Legendre, Green, HybVal>>(jParams, data);
        else
            throw std::runtime_error("addObservablesA: unknown green basis option");
        
        if( ! (jParams.is("density matrix precise") ? jParams("density matrix precise").boolean() : false))
            obs.add<DensityMatrix<Alloc, HybVal, DMStatic>>(jParams, data);
        
        if(data.dyn() != nullptr) obs.add<DensityMatrix<Alloc, HybVal, DMDynamic>>(jParams, data);
        
        if(data.occ() != nullptr && data.bullaOcc() != nullptr)
            obs.add<SuscMatrix<Alloc, true, true>>(jParams, data);
        else  if(data.occ() != nullptr)
            obs.add<SuscMatrix<Alloc, true, false>>(jParams, data);
        else if(data.bullaOcc() != nullptr)
            obs.add<SuscMatrix<Alloc, false, true>>(jParams, data);
    }
    
    
    template<typename Alloc, typename HybVal>
    void addObservablesB(jsx::value const& jParams, data::Data const& data, Observables& obs) {
        using namespace obs;
        
        obs.add<SectorProb<Alloc>>(jParams, data);
        
        if((jParams.is("quantum number susceptibility") ? jParams("quantum number susceptibility").boolean() : false) || data.bullaOcc() != nullptr)
            obs.add<SuscFlavor>(jParams, data);
        
        if((jParams.is("density matrix precise") ? jParams("density matrix precise").boolean() : false) && data.bullaOps() != nullptr)
            obs.add<Chain<Alloc, HybVal, true, true>>(jParams, data);
        else if(jParams.is("density matrix precise") ? jParams("density matrix precise").boolean() : false)
            obs.add<Chain<Alloc, HybVal, true, false>>(jParams, data);
        else if(data.bullaOps() != nullptr)
            obs.add<Chain<Alloc, HybVal, false, true>>(jParams, data);
        
        if(data.bullaOps() != nullptr) {
            if(jParams.is("green basis") ? jParams("green basis").string() == "matsubara" : true) {
                obs.add<OneParticle<BinMoments, BullaL, HybVal>>(jParams, data);
                obs.add<OneParticle<BinMoments, BullaR, HybVal>>(jParams, data);
            } else if(jParams("green basis").string() == "legendre") {
                obs.add<OneParticle<Legendre, BullaL, HybVal>>(jParams, data);
                obs.add<OneParticle<Legendre, BullaR, HybVal>>(jParams, data);
            } else
                throw std::runtime_error("obs::Observables: unknown green basis option");
        }
    }
}

#endif
