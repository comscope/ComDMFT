#ifndef EVALSIM_WORM_GREEN
#define EVALSIM_WORM_GREEN

#include "../../include/JsonX.h"
#include "../../ctqmc/include/config/Worms.h"

#include "include/Green.h"
#include "include/functions/Functions.h"
#include "include/functions/Measurements.h"
#include "include/functions/Utilities.h"

namespace evalsim {
    
    namespace worm {
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::green::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
            
            
            auto const name = cfg::green::Worm::name();
            jsx::value jWorm = jParams(name);
            
            ////Initialization
            jParams["hloc"] = ga::read_hloc<Value>("hloc.json");
            jParams["operators"] = ga::construct_annihilation_operators<Value>(jParams("hloc"));
            
            double const beta = jParams("beta").real64();
            func::iOmega const iomega(beta); auto const oneBody = jsx::at<io::Matrix<Value>>(jParams("hloc")("one body"));
            jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
            
            
            
            ////Hybridization function
            
            std::vector<io::cmat> hyb; std::vector<io::Matrix<Value>> hybMoments;
            
            std::tie(hyb, hybMoments) = partition::func::get_hybridisation<Value>(jParams);
            
            
            
            std::cout << "Reading green function ... " << std::flush;
            
            std::vector<io::cmat> green = meas::read_matrix_functions<Value,Fermion>(jMeasurements, jParams, jWorm, jHybMatrix, hyb.size());
            
            std::cout << "OK" << std::endl;
            
            
            std::cout << "Calculating self-energy with dyson ... " << std::flush;
            
            std::vector<io::cmat> self = partition::func::get_self_dyson<Value>(jParams,green,hyb);
            
            std::cout << "OK" << std::endl;
            
            /*
            std::cout << "Calculating green moments ... " << std::flush;
            
            std::vector<io::Matrix<Value>> greenMoments = func::green::compute_green_moments<Value>(jParams, hybMoments, jPartition, jObservables("partition")("scalar"));

            std::cout << "OK" << std::endl;
            
            
            std::cout << "Calculating self-energy moments ... " << std::flush;
            
            std::vector<io::Matrix<Value>> selfMoments = func::green::compute_self_moments<Value>(jParams, hybMoments, greenMoments);
            
            std::cout << "OK" << std::endl;

            
            jsx::value jObservablesOut;
            
            std::cout << "Adding self-energy high frequency tail ... "  << std::flush;
            
            func::green::add_self_tail(jHybMatrix, iomega, self, selfMoments, hyb.size());
            jObservablesOut["self-energy"] =  func::write_functions(jParams, jHybMatrix, self, selfMoments);
            
            
            std::cout << "Ok" << std::endl;
            
            
            std::cout << "Adding green function high frequency tail ... " << std::flush;
            
            func::green::add_green_tail<Value>(jParams, iomega, oneBody, hyb, self, green);
            
            jObservablesOut["green"] = func::write_functions(jParams, jHybMatrix, green, greenMoments);
               */
            
            std::cout << "Writing results ... " << std::flush;
            
            
            jsx::value jObservablesOut;
            jObservablesOut["green"] = func::write_functions<Value>(jParams, jHybMatrix, green);
            jObservablesOut["self-energy"] =  func::write_functions<Value>(jParams, jHybMatrix, self);
            
            std::cout << "Ok" << std::endl;
            
            return jObservablesOut;
        }
        
    }
    
}


#endif









