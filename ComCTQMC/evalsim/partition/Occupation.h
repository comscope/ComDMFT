#ifndef EVALSIM_PARTITION_OCCUPATION_H
#define EVALSIM_PARTITION_OCCUPATION_H


#include "ReadDensityMatrix.h"

#include "../../include/linalg/Operators.h"
#include "../../include/JsonX.h"
#include "../../include/io/Vector.h"
#include "../../include/io/Matrix.h"

namespace evalsim {
    
    namespace partition {
        
        
        template<typename Value>
        jsx::value get_occupation(jsx::value const& jParams, jsx::value const& jMeasurements, io::Matrix<Value>& occupation, io::Matrix<Value>& correlation)
        {
            jsx::value jOccupation;
            
            
            mpi::cout << "Calculating occupation ... " << std::flush;
            
            jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
            jsx::value jDensityMatrix = meas::read_density_matrix<Value>(jParams, jMeasurements("density matrix"));
        
            
            jsx::value jn = jsx::array_t(jHybMatrix.size());
            
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                linalg::mult<Value>('c', 'n', 1., jParams("operators")(i), jParams("operators")(i), .0, jn(i));
            
            std::size_t const N = jHybMatrix.size();
            std::size_t const size = N*N;
            std::size_t const rank = mpi::rank();
            std::size_t const chunk = (size + mpi::number_of_workers() - 1)/mpi::number_of_workers();
            std::vector<Value> occupation_tmp(size,0);
            std::vector<Value> correlation_tmp(size,0);
            
            for(std::size_t index = rank*chunk; index < chunk*(rank + 1); ++index){
                
                if (index>=size) break;
                
                std::size_t const i = index%(N);
                std::size_t const j = index/(N);
                
                std::string entry = jHybMatrix(i)(j).string();
                    
                if(!entry.empty()) {
                    jsx::value jOccupation;
                    linalg::mult<Value>('c', 'n', 1., jParams("operators")(i), jParams("operators")(j), .0, jOccupation);
                    occupation_tmp[i*N+j] = linalg::trace<Value>(jDensityMatrix, jOccupation);
                }
                    
                jsx::value jCorrelation;
                linalg::mult<Value>('n', 'n', 1., jn(i), jn(j), .0, jCorrelation);
                correlation_tmp[i*N+j] = linalg::trace<Value>(jDensityMatrix, jCorrelation);
            }
            
            mpi::reduce<mpi::op::sum>(occupation_tmp, mpi::master);
            mpi::reduce<mpi::op::sum>(correlation_tmp, mpi::master);
            
            for (std::size_t i=0; i<N; i++)
                for (std::size_t j=0; j<N; j++){
                    occupation(i,j) = occupation_tmp[i*N+j];
                    correlation(i,j) = correlation_tmp[i*N+j];
                }
            
            auto occupation_entries = func::get_entries(occupation, jHybMatrix);
            
            for(auto const& entry : occupation_entries)
                jOccupation[entry.first] = io::Vector<Value>{{entry.second}};
            
            occupation = func::get_matrix(occupation_entries, jHybMatrix);
            
            mpi::cout << "Ok" << std::endl;
            
            
            return jOccupation;
        }
        
    }
}


#endif









