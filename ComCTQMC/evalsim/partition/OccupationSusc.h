#ifndef EVALSIM_PARTITION_OCCUPATIONSUSC_H
#define EVALSIM_PARTITION_OCCUPATIONSUSC_H


#include "ReadDensityMatrix.h"
#include "ReadHamiltonian.h"

#include "../../include/linalg/Operators.h"
#include "../../include/JsonX.h"
#include "../../include/io/Vector.h"
#include "../../include/io/Matrix.h"
#include "../../include/options/Options.h"
#include "../../include/atomic/Generate.h"

namespace evalsim {
    
    namespace partition {
        
        
        template<typename Value>
        io::rmat get_occupation_susc_moments(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements, jsx::value& jObservables)
        {
            double const beta = jParams("beta").real64();
            jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
            
            io::rmat moments(jHybMatrix.size(), jHybMatrix.size());
            
            jsx::value jDensityMatrix = meas::read_density_matrix<Value>(jParams, jMeasurements("density matrix"));
            jsx::value jHamiltonian = get_hamiltonian<Value>(jParams);
            
            jsx::value jn = jsx::array_t(jHybMatrix.size()), jC = jsx::array_t(jHybMatrix.size());
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i) {
                linalg::mult<Value>('c', 'n', 1., jParams("operators")(i), jParams("operators")(i), .0, jn(i));
                
                linalg::mult<Value>('n', 'n',  1., jHamiltonian, jn(i), .0, jC(i));
                linalg::mult<Value>('n', 'n', -1., jn(i), jHamiltonian, 1., jC(i));
            }
            
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                    jsx::value moment;
                    
                    linalg::mult<Value>('n', 'n',  1., jC(i), jn(j), .0, moment);
                    linalg::mult<Value>('n', 'n', -1., jn(j), jC(i), 1., moment);
                    
                    moments(i, j) = ut::real(linalg::trace<Value>(jDensityMatrix, moment));
                    if(i == j) moments(i, j) -= (jsx::at<io::rvec>(jMeasurements("flavor k"))[2*i] + jsx::at<io::rvec>(jMeasurements("flavor k"))[2*i + 1])/beta;
                }
            
            return moments;
        }
        
        
        template<typename Value>
        jsx::value get_occupation_susc_bulla(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements, io::rmat const& moments, io::Matrix<Value> const& occupation, io::Matrix<Value> const& correlation, jsx::value& jObservables)
        {
            double const beta = jParams("beta").real64();
            jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");

            
            std::vector<std::vector<io::rvec>> susceptibilities(jHybMatrix.size(), std::vector<io::rvec>(jHybMatrix.size()));
            
            auto constants = moments;
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                    io::rvec measA = jsx::at<io::rvec>(jMeasurements("susceptibility flavor")(std::to_string(i) + "_" + std::to_string(j)));
                    io::rvec measB = jsx::at<io::rvec>(jMeasurements("susceptibility bulla")(std::to_string(i) + "_" + std::to_string(j)));
                    
                    io::rvec function(measA.size(), .0);
                    
                    function[0] = beta*ut::real(correlation(i, j) - occupation(i, i)*occupation(j, j)) + (measB[0] + measA[0]/beta)/2.;
                    
                    if(i == j) constants(i, j) += (jsx::at<io::rvec>(jMeasurements("flavor k"))[2*i] + jsx::at<io::rvec>(jMeasurements("flavor k"))[2*i + 1])/beta;
                    
                    for(std::size_t n = 1; n < function.size(); ++n)
                        function[n] = -beta*beta*(constants(i, j) - measB[n] - measA[n]/beta)/(4*M_PI*M_PI*n*n);
                    
                    susceptibilities[i][j] = std::move(function);
                }
            
            
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                for(std::size_t j = 0; j < jHybMatrix.size(); ++j){
                    double const moment = moments(i, j);
                    auto& function = susceptibilities[i][j];
                    std::size_t const nFit = std::max(static_cast<int>(function.size()/10.), 1);
                    
                    double A = .0, B = .0;
                    for(std::size_t n = function.size() - nFit; n < function.size(); ++n) {
                        A += moment*function[n] + (4*M_PI*M_PI*n*n/(beta*beta))*function[n]*function[n]; B += function[n]*function[n];
                    }
                    double const alpha = -A/B;
                    
                    std::size_t const nTail = std::max(static_cast<int>(beta*jPartition("susceptibility tail").real64()/(2*M_PI)), 1);
                    for(std::size_t n = function.size(); n < nTail; ++n)
                        function.push_back(-moment/((4*M_PI*M_PI*n*n/(beta*beta)) + alpha));
                }
            
            
            jsx::value jSusc;
            
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                    jSusc[std::to_string(i) + "_" + std::to_string(j)]["function"] = susceptibilities[i][j];
                    jSusc[std::to_string(i) + "_" + std::to_string(j)]["moment"]   = io::rvec{{ moments(i, j) }};
                }
            
            return jSusc;
        }
        
        
        template<typename Value>
        jsx::value get_occupation_susc_direct(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements, io::rmat const& moments, io::Matrix<Value> const& occupation, io::Matrix<Value> const& correlation, jsx::value& jObservables)
        {
            double const beta = jParams("beta").real64();
            jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
            
            std::vector<std::vector<io::rvec>> susceptibilities(jHybMatrix.size(), std::vector<io::rvec>(jHybMatrix.size()));
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                    susceptibilities[i][j] = jsx::at<io::rvec>(jMeasurements("susceptibility direct")(std::to_string(i) + "_" + std::to_string(j)));
                    susceptibilities[i][j][0] -= ut::real(occupation(i, i)*occupation(j, j))*beta;
                }
            
            jsx::value jSusc;
            
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                    jSusc[std::to_string(i) + "_" + std::to_string(j)]["function"] = susceptibilities[i][j];
                    jSusc[std::to_string(i) + "_" + std::to_string(j)]["moment"]   = io::rvec{{ moments(i, j) }};
                }
            
            return jSusc;
        }
        
    }
    
}


#endif









