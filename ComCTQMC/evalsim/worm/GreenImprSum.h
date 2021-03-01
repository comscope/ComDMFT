#ifndef EVALSIM_WORM_GREEN_IMPRSUM
#define EVALSIM_WORM_GREEN_IMPRSUM

#include "../../include/JsonX.h"
#include "../../ctqmc/include/config/Worms.h"

#include "include/functions/Functions.h"
#include "include/functions/Measurements.h"
#include "include/functions/Utilities.h"

namespace evalsim {
    
    namespace worm {
        
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::green_imprsum::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
            
            auto const name = cfg::green_imprsum::Worm::name();
            jsx::value jWorm = jParams(name);
            
            
            ////Initialization
            
            jParams["hloc"] = ga::read_hloc<Value>("hloc.json");
            
            jParams["operators"] = ga::construct_annihilation_operators<Value>(jParams("hloc"));
            
            double const beta = jParams("beta").real64();
            func::iOmega const iomega(beta); auto const oneBody = jsx::at<io::Matrix<Value>>(jParams("hloc")("one body"));
            jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
            
            std::vector<io::cmat> hyb; std::vector<io::Matrix<Value>> hybMoments;
            
            std::tie(hyb, hybMoments) = partition::func::get_hybridisation<Value>(jParams);
            
            
            std::cout << "Reading greensigma function and computing green's function... " << std::flush;
            
            std::vector<io::cmat> sigmagreen = meas::read_matrix_functions<Value,Fermion>(jMeasurements, jParams, jWorm, jHybMatrix, hyb.size());
            
            std::vector<io::cmat> green(sigmagreen.size(), io::cmat(jHybMatrix.size(), jHybMatrix.size()));
            func::green::compute_green_from_improved<Value>(jParams, iomega, oneBody, hyb, sigmagreen, green);
            
            std::cout << "OK" << std::endl;
            
            
            std::cout << "Calculating self-energy with dyson ... " << std::flush;
            
            std::vector<io::cmat> self(green.size(), io::cmat(jHybMatrix.size(), jHybMatrix.size()));
            func::green::compute_self_from_improved<Value>(jParams,sigmagreen,green,self);
            
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
            
            std::cout << "Ok" << std::endl;
             */
            
            
            
            jsx::value jObservablesOut;
            jObservablesOut["self-energy"] =  func::write_functions<Value>(jParams, jHybMatrix, self);
            jObservablesOut["green"] =  func::write_functions<Value>(jParams, jHybMatrix, green);
            
            
            return jObservablesOut;
        }
        
    }
    
}


#endif









