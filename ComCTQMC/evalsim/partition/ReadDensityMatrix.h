#ifndef EVALSIM_PARTITION_READDENSITYMATRIX_H
#define EVALSIM_PARTITION_READDENSITYMATRIX_H

#include <vector>
#include <complex>
#include <fstream>
#include <valarray>
#include <cstring>
#include <random>


#include "../../include/JsonX.h"
#include "../../include/io/Vector.h"
#include "../../include/io/Matrix.h"


namespace evalsim {
    
    namespace partition {
        
        namespace meas {

            template<typename Value> // one coud get rid of this fuck when properly implementing meas::Matrix<Value> ....
            inline jsx::value read_matrix(io::Vector<Value> const& source, int const dim) {
                if(dim*dim != source.size())
                    throw std::runtime_error("read_matrix: missmatch in size");
                
                io::Matrix<Value> dest(dim, dim);
                for(std::size_t i = 0; i < dim; ++i)
                    for(std::size_t j = 0; j < dim; ++j)
                        dest(i, j) = source.at(i*dim + j);
                
                return dest;
                
            };
            
            template<typename Value>
            inline jsx::value read_density_matrix(jsx::value const& jParams, jsx::value const& jDensityMatrixMeas) {
                jsx::value const& jEigenValues = jParams("hloc")("eigen values");
                
                jsx::value jDensityMatrix = jsx::array_t(jEigenValues.size());
                
                for(std::size_t sec = 0; sec < jEigenValues.size(); ++sec) {
                    jDensityMatrix(sec)["target"] = jsx::int64_t(sec);
                    jDensityMatrix(sec)["matrix"] = read_matrix(jsx::at<io::Vector<Value>>(jDensityMatrixMeas(sec)), jsx::at<io::rvec>(jEigenValues(sec)).size());
                }

                return jDensityMatrix;
            };
            
        }
        
    }
    
}

#endif
