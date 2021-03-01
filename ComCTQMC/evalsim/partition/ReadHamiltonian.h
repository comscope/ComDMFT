#ifndef EVALSIM_PARTITION_READHAMILTONIAN_H
#define EVALSIM_PARTITION_READHAMILTONIAN_H

#include <complex>
#include <vector>
#include <map>


#include "../../include/linalg/Operators.h"
#include "../../include/JsonX.h"
#include "../../include/io/Vector.h"
#include "../../include/io/Matrix.h"

namespace evalsim {
    
    namespace partition {
        
        template<typename Value>
        jsx::value get_hamiltonian(jsx::value const& jParams) {
            auto mu = jParams("mu").real64();
            auto filling = jsx::at<io::rvec>(jParams("hloc")("filling"));
            auto jHamiltonian = jParams("hloc")("eigen values");
            
            for(std::size_t sector = 0; sector < jHamiltonian.size(); ++sector)
                for(auto& energy : jsx::at<io::rvec>(jHamiltonian(sector)))
                    energy += -mu*filling[sector];

            return linalg::diag_to_operator<Value>(jHamiltonian);
        };
        
        
        template<typename Value>
        jsx::value get_effective_hamiltonian(jsx::value const& jParams) {
            auto mu = jParams("mu").real64();
            auto filling = jsx::at<io::rvec>(jParams("hloc")("filling"));
            auto jEffectiveHamiltonian = jParams("hloc")("eigen values");
            
            std::vector<double> shift(jEffectiveHamiltonian.size(), .0);
            
            if(jParams.is("dyn")) {
                jsx::value jq = jParams("dyn")("quantum numbers");
                jsx::value jDynFunctions = jParams("dyn")("functions");
                jsx::value jDynMatrix = jParams("dyn")("matrix");
                
                std::vector<std::vector<double>> Q(jq.size());
                std::vector<std::vector<double>> D0(jq.size(), std::vector<double>(jq.size()));
                
                if(jDynMatrix.size() != jq.size())
                    throw std::runtime_error("imp::Simple: matrix has wrong size !");
                
                for(int I = 0; I < jq.size(); ++I)
                {
                    Q[I] = ga::construct_sector_qn(jParams("hloc"), jq(I));
                    
                    if(jDynMatrix(I).size() != jq.size())
                        throw std::runtime_error("imp::Simple: matrix has wrong size !");
                    
                    for(int J = 0; J < jq.size(); ++J)
                        D0[I][J] = jsx::at<io::rvec>(jDynFunctions(jDynMatrix(I)(J).string()))[0];
                }

                for(std::size_t sector = 0; sector < jEffectiveHamiltonian.size(); ++sector)
                    for(int I = 0; I < jq.size(); ++I)
                        for(int J = 0; J < jq.size(); ++J)
                            shift[sector] += .5*Q[I][sector]*D0[I][J]*Q[J][sector];
            }

            for(std::size_t sector = 0; sector < jEffectiveHamiltonian.size(); ++sector)
                for(auto& energy : jsx::at<io::rvec>(jEffectiveHamiltonian(sector)))
                    energy += -mu*filling[sector] + shift[sector];

            return linalg::diag_to_operator<Value>(jEffectiveHamiltonian);
        };
        
    }
}


#endif









