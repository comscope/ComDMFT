#ifndef CTQMC_INCLUDE_PARAMS_H
#define CTQMC_INCLUDE_PARAMS_H

#include "Utilities.h"
#include "config/Worms.h"

#include "../../include/mpi/Utilities.h"
#include "../../include/atomic/Generate.h"
#include "../../include/options/Options.h"

namespace params {

    template<typename Value>
    void complete_impurity(jsx::value& jParams)
    {
        opt::complete_hloc<Value>(jParams);
        jParams("hloc") = ga::construct_hloc<Value>(jParams("hloc"));
        mpi::write(jParams("hloc"), "hloc.json");
        
        jParams["operators"] = ga::construct_annihilation_operators<Value>(jParams("hloc"));
        
        jParams("hybridisation")("functions") = mpi::read(jParams("hybridisation")("functions").string());
        
        if(jParams.is("dyn"))
            jParams("dyn")("functions") = mpi::read(jParams("dyn")("functions").string());
    };
    
    
    void complete_worm(jsx::value& jParams, std::string const worm)
    {
        if(jParams.is(worm)) {
            jsx::value jEntry = jParams(worm);
            jParams.object().erase(worm);
            
            jsx::value jMeas = jEntry("meas");
            jEntry.object().erase("meas");
            
            std::map<std::string, std::string> map{{"",""}, {"impr", " impr"}, {"imprsum", " imprsum"}};
            
            for(auto& meas : jMeas.array()) {
                if(!map.count(meas.string()))
                    throw std::runtime_error("params::complete_worms: invalid meas option " + meas.string());
                meas = map.at(meas.string());
            }
            
            for(auto meas : jMeas.array()) {
                jParams[worm + meas.string()] = jEntry;
                if(meas.string() != "" && jParams.is("dyn"))
                    jParams[worm] = jEntry;
            }
            
            for(auto meas : jMeas.array()) {
                jParams(worm + meas.string())["static"] = jsx::null_t();
                if(meas.string() != "" && jParams.is("dyn"))
                    jParams(worm)["dynamic"] = jsx::null_t();
            }
        }
    }
    
    void complete_worms(jsx::value& jParams)
    {
        complete_worm(jParams, cfg::green::name);
        
        complete_worm(jParams, cfg::vertex::name);
        
        complete_worm(jParams, cfg::hedin_ph::name);
        
        complete_worm(jParams, cfg::hedin_pp::name);
        
        if(jParams.is(cfg::susc_ph::name))
           jParams(cfg::susc_ph::name)["static"] = jsx::null_t();
        
        if(jParams.is(cfg::susc_pp::name))
            jParams(cfg::susc_pp::name)["static"] = jsx::null_t();
    };
}

#endif
