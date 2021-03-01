#ifndef EVALSIM_EVALSIM_H
#define EVALSIM_EVALSIM_H


#include "partition/Evalsim.h"
#include "worm/Evalsim.h"

#include "../ctqmc/include/config/Worms.h"
#include "../ctqmc/include/Params.h"
#include "../include/JsonX.h"
#include "../include/measurements/Measurements.h"
#include "../include/measurements/Error.h"
#include "../include/io/Vector.h"
#include "../include/io/Tag.h"
#include "../include/mpi/Utilities.h"


namespace evalsim {
    
    namespace worm {
        
        // catch worm-evalsims that are not yet implemented
        template<typename Value, typename W>
        jsx::value evalsim(ut::wrap<W>, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
            throw std::runtime_error("evalsim for " + W::name() + " worm not implemented");
        }
        
    }
    
    
    template<typename Value>
    bool add_dynamic(jsx::value& jStatic, jsx::value const& jDynamic)
    {
        if(jStatic.is<io::Vector<Value>>() && jDynamic.is<io::Vector<Value>>()) {
            auto& dyn = jsx::at<io::Vector<Value>>(jDynamic);
            auto& stc = jsx::at<io::Vector<Value>>(jStatic);
            
            if(stc.size() != dyn.size()) return false;
            
            for(std::size_t i = 0; i < dyn.size(); ++i) stc[i] += dyn[i];
            
            return true;
        }
        return false;
    }

    void add_dynamics(jsx::value const& jParams, jsx::value& jMeasurements, std::string const worm, std::string const meas)
    {
        if(jParams.is("dyn") && jParams.is(worm + meas)) {
            auto& jDynamic = jMeasurements(worm)("dynamic");
            auto& jStatic = jMeasurements(worm + meas)("static");

            std::set<std::string> entries;
        
            for(auto entry : jDynamic.object()) entries.insert(entry.first);
            for(auto entry : jStatic.object())  entries.insert(entry.first);
        
            for(auto entry : entries) {
                if(jDynamic.is(entry) && !jStatic.is(entry))
                    jStatic[entry] = jDynamic(entry);
                else if (jDynamic.is(entry) && jStatic.is(entry)) {
                    if(!add_dynamic<double>(jStatic(entry), jDynamic(entry)) && !add_dynamic<ut::complex>(jStatic(entry), jDynamic(entry)))
                        throw std::runtime_error("evalsim::add_dynamics: missmatch for " + worm + meas);
                }
            }
        }
    }

    
    template<typename Value>
    struct worm_clean_functor {
        template<typename W>
        void operator()(ut::wrap<W> w, jsx::value const& jParams, jsx::value& jMeasurements) const {
            if( W::name() != cfg::partition::Worm::name() && jParams.is(W::name()) ) {
                jsx::value temp = std::move(jMeasurements(W::name())("static"));
                
                jMeasurements(W::name()) = std::move(temp);
            }
        }
    };
    
    
    template<typename Value>
    struct worm_evalsim_functor {
        template<typename W>
        void operator()(ut::wrap<W> w, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value& jObservables) const {
            if( W::name() != cfg::partition::Worm::name() && jParams.is(W::name()) ) {
                std::cout << std::endl << "Begin evaluating " + W::name() + " worm measurements" << std::endl;
                
                jObservables[W::name()] = worm::evalsim<Value>(w, jParams, jMeasurements(W::name()), jMeasurements(cfg::partition::Worm::name()), jObservables);
                
                std::cout << "End evaluating " + W::name() + " worm measurements" << std::endl;
            }
        }
        void operator()(ut::wrap<cfg::partition::Worm> w, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jObservables) const {
        };
    };
    
    
    template<typename Value>
    jsx::value evalsim(jsx::value const& jParams, jsx::value& jMeasurements) {
        jsx::value jObservables;
        
        
        mpi::cout << "Begin evaluating partition measurements" << std::endl;
        
        jObservables[cfg::partition::Worm::name()] = partition::evalsim<Value>(jParams, jMeasurements(cfg::partition::Worm::name()));
        
        mpi::cout << "End evaluating partition measurements" << std::endl;
        
        add_dynamics(jParams, jMeasurements, cfg::green::name, " impr");
        add_dynamics(jParams, jMeasurements, cfg::green::name, " imprsum");
            
        add_dynamics(jParams, jMeasurements, cfg::vertex::name, " impr");
        add_dynamics(jParams, jMeasurements, cfg::vertex::name, " imprsum");
            
        add_dynamics(jParams, jMeasurements, cfg::hedin_ph::name, " impr");
        add_dynamics(jParams, jMeasurements, cfg::hedin_ph::name, " imprsum");
            
        add_dynamics(jParams, jMeasurements, cfg::hedin_pp::name, " impr");
        add_dynamics(jParams, jMeasurements, cfg::hedin_pp::name, " imprsum");
            
        cfg::for_each_type<cfg::Worm>::apply(worm_clean_functor<Value>(), jParams, jMeasurements);
            
        cfg::for_each_type<cfg::Worm>::apply(worm_evalsim_functor<Value>(), jParams, jMeasurements, jObservables);
            
        if (jParams.is("kernels")){
            jObservables["kernels"] = worm::evaluateKernels<Value>(jParams,jObservables);
            jObservables["Asymptotic Full Vertex"] = worm::evaluateFullVertexFromKernels<Value>(jParams,jObservables);
        }
        
        
        return jObservables;
        
    }
    
}

#endif












