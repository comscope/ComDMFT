#ifndef EVALSIM_PARTITION_READFUNCTIONS_H
#define EVALSIM_PARTITION_READFUNCTIONS_H

#include <vector>
#include <complex>
#include <fstream>
#include <valarray>
#include <cstring>
#include <random>

#include "Bessel.h"
#include "Functions.h"
#include "../../include/JsonX.h"
#include "../../include/io/Vector.h"
#include "../../include/io/Matrix.h"


namespace evalsim {
    
    namespace partition {
        
        namespace meas {
            
            //---------------------------------------------------------------------------------------------------------------------------------------------------------------
            //------------------------------------------------------- Functions -------------------------------------------------------------------------------------------
            //---------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            namespace impl {
                
                template<typename Value>
                inline std::pair<io::cvec, io::cvec> read_matsubara_basis(io::Vector<Value> const& data, jsx::value const& jParams, jsx::value const& jPartition, int hybSize) {
                    std::size_t const nMat = std::min(std::max(static_cast<int>(jParams("beta").real64()*jPartition("green matsubara cutoff").real64()/(2*M_PI)), 1), hybSize);
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
                
                
                template<typename Value>
                inline std::pair<io::cvec, io::cvec> read_legendre_basis(io::Vector<Value> const& data, jsx::value const& jParams, jsx::value const& jPartition, int hybSize) {
                    std::size_t const nMat = jPartition.is("green matsubara cutoff") ? std::min(std::max(static_cast<int>(jParams("beta").real64()*jPartition("green matsubara cutoff").real64()/(2*M_PI)), 1), hybSize) : hybSize;
                    std::size_t const nPol = std::min(static_cast<std::size_t>(jPartition("green legendre cutoff").int64()), data.size());
                    
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
                
                
                template<typename Value>
                inline std::pair<io::cvec, io::cvec> read_function(jsx::value const& jFunction, jsx::value const& jParams, jsx::value const& jPartition, int hybSize) {
                    if(jPartition.is("green basis") ? jPartition("green basis").string() == "matsubara" : true) {
                        return read_matsubara_basis(jsx::at<io::Vector<Value>>(jFunction), jParams, jPartition, hybSize);
                    } else if(jPartition("green basis").string() == "legendre") {
                        return read_legendre_basis(jsx::at<io::Vector<Value>>(jFunction), jParams, jPartition, hybSize);
                    } else
                        throw std::runtime_error("read_function: unknown green basis option");
                    
                    throw std::runtime_error("read_function: invalid format");
                };
                
                
                template<typename Value>
                std::pair<std::vector<io::cmat>, std::vector<io::cmat>> read_functions(jsx::value const& jFunctions, jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jHybMatrix, int hybSize)
                {
                    std::map<std::string, io::cvec> functionMatrixPos, functionMatrixNeg;
                    
                    for(auto const& function : jFunctions.object()) {
                        auto temp = impl::read_function<Value>(function.second, jParams, jPartition, hybSize);
                        
                        functionMatrixPos[function.first] = temp.first;
                        functionMatrixNeg[function.first] = temp.second;
                    }
                    
                    return std::make_pair(func::get_function_matrix(functionMatrixPos, jHybMatrix), func::get_function_matrix(functionMatrixNeg, jHybMatrix));
                };
                
                
                std::vector<io::cmat> symmetrize_functions(std::vector<io::cmat> const& functions, std::vector<io::cmat> const& functionsConj, jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jHybMatrix, int hybSize)
                {
                    std::vector<io::cmat> symmetrized(functions.size(), io::cmat(jHybMatrix.size(), jHybMatrix.size()));
                    
                    for(std::size_t n = 0; n < symmetrized.size(); ++n)
                        for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                            for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                                symmetrized[n](i, j) = (functions[n](i, j) + std::conj(functionsConj[n](j, i)))/2.;

                    return symmetrized;
                };
                
            }
            
            
            template<typename Value>
            std::vector<io::cmat> read_functions(jsx::value const& jFunctions, jsx::value const& jParams, jsx::value const& jPartition, int hybSize)
            {
                jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
                
                auto functions = impl::read_functions<Value>(jFunctions, jParams, jPartition, jHybMatrix, hybSize);

                return impl::symmetrize_functions(functions.first, functions.second, jParams, jPartition, jHybMatrix, hybSize);
            };
            
            
            template<typename Value>
            std::pair<std::vector<io::cmat>, std::vector<io::cmat>> read_conj_functions(jsx::value const& jFunctions, jsx::value const& jFunctionsConj, jsx::value const& jParams, jsx::value const& jPartition, int hybSize)
            {
                jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
                
                auto functions = impl::read_functions<Value>(jFunctions, jParams, jPartition, jHybMatrix, hybSize);
                auto functionsConj = impl::read_functions<Value>(jFunctionsConj, jParams, jPartition, jHybMatrix, hybSize);
                
                auto symmetrized = impl::symmetrize_functions(functions.first, functionsConj.second, jParams, jPartition, jHybMatrix, hybSize);
                auto symmetrizedConj = impl::symmetrize_functions(functionsConj.first, functions.second, jParams, jPartition, jHybMatrix, hybSize);
                
                return std::make_pair(symmetrized, symmetrizedConj);
            };
            
        }
        
    }
    
}

#endif
