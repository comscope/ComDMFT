#ifndef READFUNCTION
#define READFUNCTION

#include <vector>
#include <complex>
#include <fstream>
#include <valarray>
#include <cstring>
#include <random>

#include "Bessel.h"
#include "../JsonX.h"
#include "../io/Vector.h"

namespace meas {

    template<typename HybVal>
    inline std::pair<io::cvec, io::cvec> read_green_basis(io::Vector<HybVal> const& data, jsx::value const& jParams, int hybSize) {
        std::size_t const nMat = std::min(std::max(static_cast<int>(jParams("beta").real64()*jParams("green matsubara cutoff").real64()/(2*M_PI)), 1), hybSize);
        std::size_t const nIt = data.size()/4;
        double const Dtau = jParams("beta").real64()/static_cast<double>(nIt);

        io::cvec functionPos;
        io::cvec functionNeg;
        
        for(std::size_t m = 0; m < nMat; ++m) {
            double omega = M_PI*static_cast<double>(2*m + 1)/jParams("beta").real64();
            double lambda = -2.*std::sin(omega*Dtau/2.)/((Dtau*omega*(1. - omega*omega*Dtau*Dtau/24.))*jParams("beta").real64());  //missing -1/beta factor
            
            std::complex<double> iomega(.0, omega);
            std::complex<double> exp(std::exp(iomega*Dtau/2.));
            std::complex<double> fact(std::exp(iomega*Dtau));
            
            std::complex<double> tempPos = .0;
            std::complex<double> tempNeg = .0;
            auto const* itData = data.data();
            for(std::size_t i = 0; i < nIt; ++i) {
                std::complex<double> coeff = lambda*exp;
                
                tempPos += coeff**itData;
                tempNeg += std::conj(coeff)**itData;
                ++itData; coeff *= iomega;
                
                tempPos += coeff**itData;
                tempNeg += std::conj(coeff)**itData;
                ++itData; coeff *= iomega/2.;
                
                tempPos += coeff**itData;
                tempNeg += std::conj(coeff)**itData;
                ++itData; coeff *= iomega/3.;
                
                tempPos += coeff**itData;
                tempNeg += std::conj(coeff)**itData;
                ++itData;
                
                exp *= fact;
            }
            
            functionPos.push_back(tempPos);
            functionNeg.push_back(tempNeg);
        }
 
        return std::make_pair(functionPos, functionNeg);
    };
    
    template<typename HybVal>
    inline std::pair<io::cvec, io::cvec> read_legendre_basis(io::Vector<HybVal> const& data, jsx::value const& jParams, int hybSize) {
        std::size_t const nMat = jParams.is("green matsubara cutoff") ? std::min(std::max(static_cast<int>(jParams("beta").real64()*jParams("green matsubara cutoff").real64()/(2*M_PI)), 1), hybSize) : hybSize;
        std::size_t const nPol = std::min(static_cast<std::size_t>(jParams("green legendre cutoff").int64()), data.size());
        
        io::cvec functionPos;
        io::cvec functionNeg;
        
        be::Bessel bessel(nPol, nMat);
        for(std::size_t m = 0; m < nMat; ++m) {
            std::complex<double> tempPos = .0;
            std::complex<double> tempNeg = .0;
            
            std::complex<double> const I(.0, 1.);
            std::complex<double> fact = m%2 ? -1. : 1.;
            
            for(std::size_t n = 0; n < nPol; ++n) {
                fact *= I;
                tempPos += fact*bessel(n, m)*data[n];
                tempNeg += std::conj(fact*bessel(n, m))*data[n];
            };
            
            functionPos.push_back(tempPos);
            functionNeg.push_back(tempNeg);
        };
        
        return std::make_pair(functionPos, functionNeg);
    };
    
    
    inline std::pair<io::cvec, io::cvec> read_function(jsx::value const& jFunction, jsx::value const& jParams, int hybSize) {
        
        if(jParams.is("green basis") ? jParams("green basis").string() == "matsubara" : true) {
            
            if(jFunction.is<io::rvec>())
                return read_green_basis(jsx::at<io::rvec>(jFunction), jParams, hybSize);
            else if(jFunction.is<io::cvec>())
                return read_green_basis(jsx::at<io::cvec>(jFunction), jParams, hybSize);
                
        } else if(jParams("green basis").string() == "legendre") {
            
            if(jFunction.is<io::rvec>())
                return read_legendre_basis(jsx::at<io::rvec>(jFunction), jParams, hybSize);
            else if(jFunction.is<io::cvec>())
                return read_legendre_basis(jsx::at<io::cvec>(jFunction), jParams, hybSize);
            
        } else
            throw std::runtime_error("meas::read_function: unknown green basis option");
        
        throw std::runtime_error("meas::read_function: invalid format");
        
    };
}

#endif
