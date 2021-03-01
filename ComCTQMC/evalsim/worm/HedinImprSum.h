#ifndef EVALSIM_WORM_HEDIN_IMPRSUM
#define EVALSIM_WORM_HEDIN_IMPRSUM

#include "../../include/JsonX.h"
#include "../../ctqmc/include/config/Worms.h"

#include "include/Common.h"
#include "include/Hedin.h"
#include "include/functions/Functions.h"
#include "include/functions/Measurements.h"
#include "include/functions/Utilities.h"

namespace evalsim {
    
    namespace worm {
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::hedin_ph_imprsum::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
            
            auto const name = cfg::hedin_ph_imprsum::Worm::name();
            jsx::value jWorm = jParams(name);
            
            ////Initialization
            func::Frequencies<Value> frequencies(jWorm);
            
            std::cout << "Reading hybridisation function ... " << std::flush;
            
            jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
            std::vector<io::cmat> hyb; std::vector<io::Matrix<Value>> hybMoments;
            
            std::tie(hyb, hybMoments) = partition::func::get_hybridisation<Value>(jParams);
            
            std::cout << "Ok" << std::endl;
            
            
            std::cout << "Reading " << name << " function ... " << std::flush;
            
            std::vector<io::ctens> hedin_impr = meas::read_tensor_functions<Value,Fermion,Boson>(jMeasurements, jParams, jWorm, jHybMatrix, hyb.size());
            
            std::cout << "Ok" << std::endl;
            
            
            std::cout << "Gathering one-particle green's and self-energy functions ... " << std::flush;
                        
            auto const green_half_axis = func::get_green_from_obs<Value>(jParams,jObservables,jHybMatrix,hyb.size(),"green");
            auto const green = func::green_function_on_full_axis( green_half_axis );
            auto const self_half_axis = func::get_green_from_obs<Value>(jParams,jObservables,jHybMatrix,hyb.size(),"self-energy");
            auto const self = func::green_function_on_full_axis( self_half_axis );
            func::OmegaMap greenOM(green.size(),false,true);
            
            std::cout << "Ok" << std::endl;
            
            
            std::cout << "Calculating susceptibility ... " << std::flush;
            
            std::vector<io::ctens> susc(hedin_impr.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            func::hedin::ph::compute_connected_from_improved_estimator<Value>(jParams, frequencies, green, self, greenOM, jObservables("partition")("occupation"), hedin_impr, susc);
            
            std::cout << "OK" << std::endl;
            
            
            std::cout << "Enforcing symmetries ... " << std::flush;
            
            std::vector<io::ctens> susc_symm(susc.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            func::hedin::ph::enforce_symmetries<Value>(jParams, frequencies, susc, susc_symm);
            
            std::cout << "Ok" << std::endl;
            
            
            jsx::value jObservablesOut;
            
            jObservablesOut["susceptibility"] = func::write_functions<Value>(jParams, jHybMatrix, susc_symm);
            
            return jObservablesOut;
        }
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::hedin_pp_imprsum::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
            
            auto const name = cfg::hedin_pp_imprsum::Worm::name();
            jsx::value jWorm = jParams(name);
            
            ////Initialization
            func::Frequencies<Value> frequencies(jWorm);
            
            std::cout << "Reading hybridisation function ... " << std::flush;
            
            jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
            std::vector<io::cmat> hyb; std::vector<io::Matrix<Value>> hybMoments;
            
            std::tie(hyb, hybMoments) = partition::func::get_hybridisation<Value>(jParams);
            
            std::cout << "Ok" << std::endl;
            
            
            std::cout << "Reading " << name << " function ... " << std::flush;
            
            std::vector<io::ctens> hedin_impr = meas::read_tensor_functions<Value,Fermion,Boson>(jMeasurements, jParams, jWorm, jHybMatrix, hyb.size());
            
            std::cout << "Ok" << std::endl;
            
            
            std::cout << "Gathering one-particle green's and self-energy functions ... " << std::flush;
            
            auto const green_half_axis = func::get_green_from_obs<Value>(jParams,jObservables,jHybMatrix,hyb.size(),"green");
            auto const green = func::green_function_on_full_axis( green_half_axis );
            auto const self_half_axis = func::get_green_from_obs<Value>(jParams,jObservables,jHybMatrix,hyb.size(),"self-energy");
            auto const self = func::green_function_on_full_axis( self_half_axis );
            func::OmegaMap greenOM(green.size(),false,true);
            
            std::cout << "Ok" << std::endl;
            
            
            std::cout << "Calculating susceptibility ... " << std::flush;
            
            std::vector<io::ctens> susc(hedin_impr.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            func::hedin::pp::compute_connected_from_improved_estimator<Value>(jParams, frequencies, green, self, greenOM, jObservables("partition")("occupation"), hedin_impr, susc);
            
            std::cout << "OK" << std::endl;
            
            
            std::cout << "Enforcing symmetries ... " << std::flush;
            
            std::vector<io::ctens> susc_symm(susc.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            func::hedin::pp::enforce_symmetries<Value>(jParams, frequencies, susc, susc_symm);
            
            std::cout << "Ok" << std::endl;
            
            
            jsx::value jObservablesOut;
            
            jObservablesOut["susceptibility"] = func::write_functions<Value>(jParams, jHybMatrix, susc_symm);
            
            return jObservablesOut;
        }
        
        
    }
    
}




#endif









