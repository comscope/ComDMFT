#ifndef READDENSITYMATRIX
#define READDENSITYMATRIX

#include <vector>
#include <complex>
#include <fstream>
#include <valarray>
#include <cstring>
#include <random>

#include "../JsonX.h"
#include "../io/Vector.h"
#include "../io/Matrix.h"
#include "../linalg/LinAlg.h"
#include "../../ctqmc/include/Utilities.h"

namespace meas {
    
    template<typename HybVal>
    inline jsx::value read_matrix(io::Vector<HybVal> const& source, int const dim) {
        if(dim*dim != source.size())
            throw std::runtime_error("meas::read_density_matrix: missmatch in size");
        
        io::cmat dest(dim, dim);
        for(std::size_t i = 0; i < dim; ++i)
            for(std::size_t j = 0; j < dim; ++j)
                dest(i, j) = source.at(i*dim + j);
        
        return dest;
        
    };
    
    inline jsx::value read_density_matrix_impl(jsx::value& jDensityMatrixMeas, jsx::value& jEigenValues) {
        jsx::value jDensityMatrix = jsx::array_t(jEigenValues.size());
        
        for(std::size_t sec = 0; sec < jEigenValues.size(); ++sec) {
            jDensityMatrix(sec)["target"] = jsx::int64_t(sec);
            
            if(jDensityMatrixMeas(sec).is<io::rvec>())
                jDensityMatrix(sec)["matrix"] = read_matrix(jsx::at<io::rvec>(jDensityMatrixMeas(sec)), jsx::at<io::rvec>(jEigenValues(sec)).size());
            else if(jDensityMatrixMeas(sec).is<io::cvec>())
                jDensityMatrix(sec)["matrix"] = read_matrix(jsx::at<io::cvec>(jDensityMatrixMeas(sec)), jsx::at<io::rvec>(jEigenValues(sec)).size());
            else
                throw std::runtime_error("meas::read_density_matrix: unknown format");
        }
        
        return jDensityMatrix;
    };
    
    
    inline jsx::value read_density_matrix(jsx::value& jParams, jsx::value& jMeasurements) {
        return read_density_matrix_impl(jMeasurements("density matrix"), jParams("hloc")("eigen values"));
    };
    
    
    inline jsx::value read_density_matrix_prime(jsx::value& jParams, jsx::value& jMeasurements, double& Q2) {
        jsx::value jDensityMatrix = read_density_matrix_impl(jMeasurements("density matrix"), jParams("hloc")("eigen values"));
        jsx::value jDensityMatrixDyn = read_density_matrix_impl(jMeasurements("density matrix dyn"), jParams("hloc")("eigen values"));
        
        double const beta = jParams("beta").real64();
        double const UrOmega0 = jParams("dyn")(0).real64(); double UrTau0 = .0;
        for(std::size_t n = 0; n < jParams("dyn").size(); ++n)
            UrTau0 += 2./beta*jParams("dyn")(n).real64();
        
        Q2 = -UrTau0 + jsx::at<io::rvec>(jMeasurements("dyndyn"))[0];
        
        jsx::value jDensityMatrixPrime = jsx::array_t(jDensityMatrix.size());
        for(std::size_t s = 0; s < jDensityMatrix.size(); ++s) {
            jDensityMatrixPrime(s)["target"] = jDensityMatrix(s)("target").real64();
            
            auto const filling = jsx::at<io::rvec>(jParams("hloc")("filling")).at(s);
            auto& block = jsx::at<io::cmat>(jDensityMatrix(s)("matrix"));
            auto& blockDyn = jsx::at<io::cmat>(jDensityMatrixDyn(s)("matrix"));
            
            Q2 -= 2.*UrOmega0*filling*linalg::trace(blockDyn).real();
            Q2 += UrOmega0*UrOmega0*filling*filling*linalg::trace(block).real();
            
            io::cmat blockPrime = blockDyn;
            for(std::size_t i = 0; i < block.I(); ++i)
                for(std::size_t j = 0; j < block.J(); ++j)
                    blockPrime(i, j) -= UrOmega0*filling*block(i, j);
            
            jDensityMatrixPrime(s)["matrix"] = blockPrime;
        }
        
        return jDensityMatrixPrime;
    };
}

#endif
