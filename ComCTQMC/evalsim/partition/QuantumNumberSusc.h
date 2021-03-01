#ifndef EVALSIM_PARTITION_QUANTUMNUMBERSUSC_H
#define EVALSIM_PARTITION_QUANTUMNUMBERSUSC_H


#include "../../include/linalg/Operators.h"
#include "../../include/JsonX.h"
#include "../../include/io/Vector.h"
#include "../../include/io/Matrix.h"
#include "../../include/options/Options.h"
#include "../../include/atomic/Generate.h"

namespace evalsim {
    
    namespace partition {
        
        
        jsx::value get_qn_susc(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements, jsx::value const& jScalar) {
            double const beta = jParams("beta").real64(); jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
            
            jsx::value jSusc;
            
            for(auto jqn : jPartition("quantum numbers").object()) {
                
                mpi::cout << "Reading " << jqn.first << " susceptibility ... " << std::flush;
                
                auto const qn = jsx::at<io::rvec>(jqn.second);
                
                double moment = 0;
                io::rvec function(jsx::at<io::rvec>(jMeasurements("susceptibility flavor")("0_0")).size(), .0);
                for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                    for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                        double const fact = qn.at(i)*qn.at(j);
                        auto const meas = jsx::at<io::rvec>(jMeasurements("susceptibility flavor")(std::to_string(i) + "_" + std::to_string(j)));
                        
                        function[0] += fact*meas[0]/(2.*beta);
                        
                        for(std::size_t n = 1; n < function.size(); ++n) function[n] += fact*beta*meas[n]/(4*M_PI*M_PI*n*n);
                        
                        if(i == j) moment -= fact*(jsx::at<io::rvec>(jMeasurements("flavor k"))[2*i] + jsx::at<io::rvec>(jMeasurements("flavor k"))[2*i + 1])/beta;
                    }
                function[0] += beta*jsx::at<io::rvec>(jScalar(jqn.first + jqn.first)).at(0);
                function[0] -= beta*jsx::at<io::rvec>(jScalar(jqn.first)).at(0)*jsx::at<io::rvec>(jScalar(jqn.first)).at(0);
                
                mpi::cout << "Ok" << std::endl;
                
                
                mpi::cout << "Adding " << jqn.first << " susceptibility high frequency tail ... " << std::flush;
                
                std::size_t const nFit = std::max(static_cast<int>(function.size()/10.), 1);
                
                double A = .0, B = .0;
                for(std::size_t n = function.size() - nFit; n < function.size(); ++n) {
                    A += moment*function[n] + (4*M_PI*M_PI*n*n/(beta*beta))*function[n]*function[n]; B += function[n]*function[n];
                }
                double const alpha = -A/B;
                
                std::size_t const nTail = std::max(static_cast<int>(beta*jPartition("susceptibility tail").real64()/(2*M_PI)), 1);
                for(std::size_t n = function.size(); n < nTail; ++n)
                    function.push_back(-moment/((4*M_PI*M_PI*n*n/(beta*beta)) + alpha));
                
                jSusc[jqn.first]["function"] = function;
                jSusc[jqn.first]["moment"] = io::rvec{{ moment }};
                
                mpi::cout << "Ok" << std::endl;
            }
            
            
            return jSusc;
        }
        
    }
}


#endif









