#ifndef EVALSIM_WORM_SUSC
#define EVALSIM_WORM_SUSC


#include "../../include/JsonX.h"
#include "../../ctqmc/include/config/Worms.h"

#include "include/Susc.h"
#include "include/functions/Functions.h"
#include "include/functions/Measurements.h"
#include "include/functions/Utilities.h"

namespace evalsim {
    
    namespace worm {
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::susc_ph::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
            
            auto const name = cfg::susc_ph::Worm::name();
            jsx::value jWorm = jParams(name);
            
            
            func::BosonFrequencies<Value> frequencies(jWorm);
            
            ////Hybridization function
            
            std::cout << "Reading hybridisation function ... " << std::flush;
            
            jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
            std::vector<io::cmat> hyb; std::vector<io::Matrix<Value>> hybMoments;
            
            std::tie(hyb, hybMoments) = partition::func::get_hybridisation<Value>(jParams);
            
            std::cout << "Ok" << std::endl;
                                  
                                  
            std::cout << "Reading " << name << " function ... " << std::flush;
            
            std::vector<io::ctens> susc = meas::read_tensor_functions<Value,Boson>(jMeasurements, jParams, jWorm, jHybMatrix, hyb.size());
            
            std::cout << "Ok" << std::endl;
                                  
            
            std::cout << "Computing susceptibility ... " << std::flush;
            
            func::susc::ph::compute_and_subtract_disconnected<Value>(jParams,jObservables("partition")("occupation"), frequencies, susc);
            
            std::cout << "Ok" << std::endl;
            
            std::cout << "Enforcing symmetries ... " << std::flush;
            
            std::vector<io::ctens> susc_symm(susc.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            func::susc::ph::enforce_symmetries<Value>(jParams,jWorm,susc,susc_symm);
            
            std::cout << "Ok" << std::endl;
            
            jsx::value jObservablesOut;
            
            jObservablesOut["susceptibility"] = func::write_functions<Value>(jParams, jHybMatrix, susc_symm);
            
            return jObservablesOut;
        }
        
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::susc_pp::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
            
            auto const name = cfg::susc_pp::Worm::name();
            jsx::value jWorm = jParams(name);
            
            ////Hybridization function
            
            std::cout << "Reading hybridisation function ... " << std::flush;
            
            jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
            std::vector<io::cmat> hyb; std::vector<io::Matrix<Value>> hybMoments;
            
            std::tie(hyb, hybMoments) = partition::func::get_hybridisation<Value>(jParams);
            
            std::cout << "Ok" << std::endl;
            
            
            std::cout << "Reading " << name << " function ... " << std::flush;
            
            std::vector<io::ctens> susc = meas::read_tensor_functions<Value,Boson>(jMeasurements, jParams, jWorm, jHybMatrix, hyb.size());
            
            std::cout << "Ok" << std::endl;
            
            // Particle particle susceptibility is equal to the two-particle green function
            std::cout << "Enforcing symmetries ... " << std::flush;
            
            std::vector<io::ctens> susc_symm(susc.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            func::susc::pp::enforce_symmetries<Value>(jParams,susc,susc_symm);
        
            std::cout << "Ok" << std::endl;
            
            jsx::value jObservablesOut;
            
            jObservablesOut["susceptibility"] = func::write_functions<Value>(jParams, jHybMatrix, susc_symm);
            
            return jObservablesOut;
        }
        
        
    }
    
}


#endif









