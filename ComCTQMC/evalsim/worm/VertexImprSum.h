#ifndef EVALSIM_WORM_VERTEX_IMPRSUM
#define EVALSIM_WORM_VERTEX_IMPRSUM


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
        jsx::value evalsim(ut::wrap<cfg::vertex_imprsum::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
            
            auto const name = cfg::vertex_imprsum::Worm::name();
            jsx::value jWorm = jParams(name);
            
            //Initialization
            func::Frequencies<Value> frequencies(jWorm);
            
            ////Hybridization function
            jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
            std::vector<io::cmat> hyb; std::vector<io::Matrix<Value>> hybMoments;
            
            std::tie(hyb, hybMoments) = partition::func::get_hybridisation<Value>(jParams);
            
            std::cout << "Reading " << name << " function ... " << std::flush;
            
            std::vector<io::ctens> greentwo_impr = meas::read_tensor_functions<Value,Fermion,Fermion,Boson>(jMeasurements, jParams, jWorm, jHybMatrix, hyb.size());
            
            std::cout << "OK" << std::endl;
            
            
            std::cout << "Gathering one-particle green's function ... " << std::flush;
            
            auto const green_half_axis = func::get_green_from_obs<Value>(jParams,jObservables,jHybMatrix,hyb.size(),"green");
            auto const green = func::green_function_on_full_axis( green_half_axis );
            auto const self_half_axis = func::get_green_from_obs<Value>(jParams,jObservables,jHybMatrix,hyb.size(),"self-energy");
            auto const self = func::green_function_on_full_axis( self_half_axis );
            func::OmegaMap greenOM(green.size(),false,true);
            
            std::cout << "Ok" << std::endl;
            
            std::cout << "Calculating connected susceptibility from improved estimator ... " << std::flush;
            
            std::vector<io::ctens> greentwo(greentwo_impr.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            
            func::vertex::compute_connected_from_improved_estimator<Value>(jParams, frequencies, green, self, greenOM, greentwo_impr, greentwo);
            
            std::cout << "Ok" << std::endl;
            
            std::cout << "Enforcing symmetries ... " << std::flush;
            
            std::vector<io::ctens> susc_symm(greentwo.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            func::vertex::enforce_symmetries<Value>(jParams, frequencies, greentwo, susc_symm);
            
            std::cout << "Ok" << std::endl;
            
            std::vector<io::ctens> full_vertex(greentwo.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            
            if(jWorm.is("full") ? jWorm("full").boolean() : false){
                
                std::cout << "Computing full vertex ... " << std::flush;
                
                func::vertex::compute_full_vertex<Value>(jParams, frequencies, green, greenOM, susc_symm, full_vertex);
                
                std::cout << "Ok" << std::endl;
            }
            
            
            std::cout << "Adding results to params.obs.json ... " << std::flush;
                
            jsx::value jObservablesOut;
            
            jObservablesOut["susceptibility"] = func::write_functions<Value>(jParams, jHybMatrix, susc_symm);
            
            if(jWorm.is("full") ? jWorm("full").boolean() : false)
                jObservablesOut["full vertex"] = func::write_functions<Value>(jParams, jHybMatrix, full_vertex);
                
            std::cout << "Ok" << std::endl;
            
            return jObservablesOut;
        }
        
    }
    
}


#endif









