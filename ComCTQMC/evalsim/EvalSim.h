#ifndef EVALSIM_H
#define EVALSIM_H

#include <complex>
#include <vector>
#include <map>

#include "Functions.h"

#include "../include/linalg/Operators.h"
#include "../include/JsonX.h"
#include "../include/io/Vector.h"
#include "../include/io/Matrix.h"
#include "../include/measurements/ReadFunction.h"

#include "../ctqmc/include/bath/Hyb.h"

struct iOmega {
    iOmega() = delete;
    iOmega(double beta) : beta_(beta) {};
    std::complex<double> operator()(int n) const {
        return {.0, (2*n + 1)*M_PI/beta_};
    };
private:
    double const beta_;
};


io::cmat get_hybridisation_moments(std::map<std::string, io::cvec> const& functions, jsx::value const& jParams, jsx::value const& jHybMatrix)
{
    std::map<std::string, ut::complex> moments;
    for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
        for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
            if(jHybMatrix(i)(j).string() != "") {
                if(jHybMatrix(j)(i).string() == "")
                    throw std::runtime_error("Hyb: invalid hybridisation matrix.");
                
                std::string const entry = jHybMatrix(i)(j).string();
                std::string const entryTransp = jHybMatrix(j)(i).string();
                
                if(!moments.count(entry)) {
                    if(jParams("complex hybridisation").boolean())
                        moments[entry] = bath::Fit<ut::complex>(jParams("beta").real64(), functions.at(entry), functions.at(entryTransp)).moment();
                    else
                        moments[entry] = bath::Fit<double>(jParams("beta").real64(), functions.at(entry), functions.at(entryTransp)).moment();
                }
            }
    
    return get_matrix(moments, jHybMatrix);
}


std::vector<io::cmat> read_functions(jsx::value& jFunctions, jsx::value& jParams, jsx::value const& jHybMatrix, int hybSize)
{
    std::map<std::string, io::cvec> functionsPos, functionsNeg;
    
    for(auto& function : jFunctions.object()) {
        auto temp = meas::read_function(function.second, jParams, hybSize);
        
        functionsPos[function.first] = temp.first;
        functionsNeg[function.first] = temp.second;
    }
    
    std::vector<io::cmat> matricesPos = get_function_matrix(functionsPos, jHybMatrix);
    std::vector<io::cmat> matricesNeg = get_function_matrix(functionsNeg, jHybMatrix);
    
    for(std::size_t n = 0; n < matricesPos.size(); ++n)
        for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
            for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                matricesPos[n](i, j) = (matricesPos[n](i, j) + std::conj(matricesNeg[n](j, i)))/2.;
    
    // here the function is symmetric if real in imaginary time
    
    return matricesPos;
};


void add_self_tail(jsx::value const& jHybMatrix, iOmega const& iomega, std::vector<io::cmat>& function, std::vector<io::cmat> const& moments, std::size_t hybSize)
{
    io::cmat alpha(jHybMatrix.size(), jHybMatrix.size());
    io::cmat gamma(jHybMatrix.size(), jHybMatrix.size());
    
    std::size_t const nFit  = std::max(function.size()/8, static_cast<std::size_t>(1));

    for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
        for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
            if(jHybMatrix(i)(j).string() != "") {
                std::complex<double> iomegaAvg = .0;
                std::complex<double> selfPosAvg = .0;
                std::complex<double> selfNegAvg = .0;
                
                for(std::size_t n = function.size() - nFit;  n < function.size(); ++n) {
                    iomegaAvg += iomega(n);
                    selfPosAvg += function[n](i, j);
                    selfNegAvg += std::conj(function[n](j, i));
                }
                
                iomegaAvg /= static_cast<double>(nFit);
                selfPosAvg /= static_cast<double>(nFit);
                selfNegAvg /= static_cast<double>(nFit);
                
                alpha(i, j) = -moments[1](i, j)/2.*(1./(selfPosAvg - moments[0](i, j)) + 1./(selfNegAvg - moments[0](i, j)));
                gamma(i, j) = iomegaAvg*moments[1](i, j)/2.*(-1./(selfPosAvg - moments[0](i, j)) + 1./(selfNegAvg - moments[0](i, j))) + iomegaAvg*iomegaAvg;
            }
    
    
    for(std::size_t n = function.size(); n < hybSize; ++n) {
        io::cmat temp(jHybMatrix.size(), jHybMatrix.size());
        
        for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
            for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                if(jHybMatrix(i)(j).string() != "")
                    temp(i, j) = moments[0](i, j) + moments[1](i, j)/(iomega(n) - alpha(i, j) - gamma(i, j)/iomega(n));
        
        function.push_back(temp);
    }
}


jsx::value write_functions(jsx::value const& jParams, jsx::value const& jHybMatrix, std::vector<io::cmat> const& functionsMatrix, std::vector<io::cmat> const& momentsMatrix) {
    jsx::value jFunction;
    
    std::map<std::string, io::cvec> functions = get_function_entries(functionsMatrix, jHybMatrix);
    for(auto& function : functions)
        jFunction[function.first]["function"] = std::move(function.second);
    
    std::map<std::string, io::cvec> moments = get_function_entries(momentsMatrix, jHybMatrix);
    
    if(jParams("complex hybridisation").boolean()) {
        for(auto& moment : moments) jFunction[moment.first]["moments"] = moment.second;
    } else {
        for(auto& moment : moments) {
            io::rvec momentReal(moment.second.size());
            for(std::size_t i = 0; i < moment.second.size(); ++i) momentReal[i] = moment.second[i].real();
            jFunction[moment.first]["moments"] = momentReal;
        }
    }
    
    return jFunction;
}


jsx::value get_hamiltonian(jsx::value& jParams) {
    auto jEigenValues = jParams("hloc")("eigen values"); auto& jN = jParams("hloc")("filling");
    
    for(std::size_t sector = 0; sector < jEigenValues.size(); ++sector)
        for(auto& energy : jsx::at<io::rvec>(jEigenValues(sector)))
            energy += -jParams("mu").real64()*jsx::at<io::rvec>(jN)[sector];
    
    jsx::value jHamiltonian = linalg::diag_to_operator(jEigenValues);
    linalg::make_operator_complex(jHamiltonian);
    return jHamiltonian;
};


jsx::value get_effective_hamiltonian(jsx::value& jParams) {
    auto jEigenValues = jParams("hloc")("eigen values"); auto& jN = jParams("hloc")("filling");
    auto const Ur0 = jParams.is("dyn") ? jParams("dyn")(0).real64() : .0;
    
    for(std::size_t sector = 0; sector < jEigenValues.size(); ++sector)
        for(auto& energy : jsx::at<io::rvec>(jEigenValues(sector)))
            energy += -jParams("mu").real64()*jsx::at<io::rvec>(jN)[sector] + .5*Ur0*jsx::at<io::rvec>(jN)[sector]*jsx::at<io::rvec>(jN)[sector];
    
    jsx::value jHamiltonian = linalg::diag_to_operator(jEigenValues);
    linalg::make_operator_complex(jHamiltonian);
    return jHamiltonian;
};


struct VecLess {
    VecLess(std::size_t size, std::size_t start, std::size_t end) : size_(size), start_(start), end_(end) {};
    bool operator()(std::vector<double> const& lhs, std::vector<double> const& rhs) const {
        if(lhs.size() != size_ || rhs.size() != size_)
            throw std::runtime_error("VecCompare");
        
        for(std::size_t i = start_; i < end_; ++i) {
            if(lhs[i] > rhs[i])
                return true;
            else if(lhs[i] < rhs[i])
                return false;
        }
        
        return false;
    }
private:
    std::size_t const size_, start_, end_;
};


double truncate(double val, int prec) {
    if(std::abs(val) < 1.e-8) return .0;
    
    std::stringstream temp;
    temp << std::setprecision(prec) << val;
    temp >> val; return val;
}




#endif









