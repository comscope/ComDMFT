#ifndef EVALSIM_WORM_VERTEX
#define EVALSIM_WORM_VERTEX

#include "../../include/JsonX.h"
#include "../../ctqmc/include/config/Worms.h"

#include "include/Common.h"
#include "include/Vertex.h"
#include "include/functions/Functions.h"
#include "include/functions/Measurements.h"
#include "include/functions/Utilities.h"

namespace evalsim {
    
    namespace worm {
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::vertex::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
            
            auto const name = cfg::vertex::Worm::name();
            jsx::value jWorm = jParams(name);
            
            //Initialization
            func::Frequencies<Value> frequencies(jWorm);
            
            ////Hybridization function
            
            jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
            std::vector<io::cmat> hyb; std::vector<io::Matrix<Value>> hybMoments;
            
            std::tie(hyb, hybMoments) = partition::func::get_hybridisation<Value>(jParams);
            
            
            std::cout << "Reading " << name << " function ... " << std::flush;
            
            std::vector<io::ctens> greentwo = meas::read_tensor_functions<Value,Fermion,Fermion,Boson>(jMeasurements, jParams, jWorm, jHybMatrix, hyb.size());
            
            std::cout << "OK" << std::endl;
            
            std::cout << "Gathering one-particle green's function ... " << std::flush;
            
            auto const green_half_axis = func::get_green_from_obs<Value>(jParams,jObservables,jHybMatrix,hyb.size(),"green");
            auto const green = func::green_function_on_full_axis( green_half_axis );
            func::OmegaMap greenOM(green.size(),false,true);
            
            std::cout << "Ok" << std::endl;
            
            
            std::cout << "Calculating connected part of two-particle green's function " << std::flush;
            
            func::vertex::compute_and_subtract_disconnected<Value>(jParams, frequencies, green, greenOM, greentwo);
            
            std::cout << "Ok" << std::endl;
            
            std::cout << "Enforcing symmetries ... " << std::flush;
            
            std::vector<io::ctens> susc_symm(greentwo.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            func::vertex::enforce_symmetries<Value>(jParams, frequencies, greentwo,susc_symm);
            
            std::cout << "Ok" << std::endl;
            
            
            std::cout << "Adding results to params.obs.json ... " << std::flush;
            
            jsx::value jObservablesOut;
            
            jObservablesOut["susceptibility"] = func::write_functions<Value>(jParams, jHybMatrix, susc_symm);
            
            std::cout << "Ok" << std::endl;
            
            return jObservablesOut;
        }
        
    }
    
}


#endif









