#ifndef EVALSIM_WORM_INCLUDE_GREEN
#define EVALSIM_WORM_INCLUDE_GREEN


#include "../../../include/JsonX.h"
#include "../../../ctqmc/include/config/Worms.h"
#include "functions/Functions.h"
#include "functions/Measurements.h"
#include "functions/Utilities.h"
#include "Common.h"


namespace evalsim {
    
    namespace worm {
        
        namespace func {
        
            namespace green {
                
                template<typename Value>
                void compute_green_from_improved(jsx::value const& jParams, iOmega const& iomega, io::Matrix<Value> const& oneBody, std::vector<io::cmat> const& hyb, std::vector<io::cmat> const& sigmagreen, std::vector<io::cmat> & green){
                    
                    double const mu = jParams("mu").real64();
                    jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
                    
                    for(std::size_t n = 0; n < sigmagreen.size(); ++n)
                        for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                            for(std::size_t j = 0; j < jHybMatrix.size(); ++j){
                                green[n](i,j)=0;
                                for(std::size_t k = 0; k < jHybMatrix.size(); ++k)
                                    green[n](i,j) += std::abs(hyb[n](i,k)) ? //
                                    ((j == k ? 1.0 : 0.0) + sigmagreen[n](j,k))/ //(delta_ik + greensigma_ik
                                    ((i == k ? iomega(n) + mu : .0) - oneBody(i, k) - hyb[n](i, k)) // G0_ik
                                    : 0 ; //DMFT assumption?
                            }
                    
                }
                
                template<typename Value>
                void compute_self_from_improved(jsx::value const& jParams, std::vector<io::cmat> const& sigmagreen, std::vector<io::cmat> const& green, std::vector<io::cmat>& self){
                    
                    jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
                    
                    for(std::size_t n = 0; n < sigmagreen.size(); ++n)
                        for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                            for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                                self[n](i,j) = std::abs(green[n](i,j)) ? sigmagreen[n](i,j)/green[n](i,j) : 0;
                    
                    
                }
                
                template <typename Value>
                std::vector<io::Matrix<Value>> compute_green_moments(jsx::value const& jParams, std::vector<io::Matrix<Value>> const& hybMoments, jsx::value const& jPartition, jsx::value const& jScalar){
                    double const beta = jParams("beta").real64();
                    jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
                    jsx::value jDensityMatrix = partition::meas::read_density_matrix<Value>(jParams, jPartition("density matrix"));
                    
                    std::vector<io::Matrix<Value>> greenMoments;
                    
                    jsx::value jHamiltonian = partition::get_hamiltonian<Value>(jParams);
                    
                    jsx::value jC = jsx::array_t(jHybMatrix.size());
                    for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                    {
                        linalg::mult<Value>('n', 'n',  1., jHamiltonian, jParams("operators")(i), .0, jC(i));
                        linalg::mult<Value>('n', 'n', -1., jParams("operators")(i), jHamiltonian, 1., jC(i));
                    }
                    
                    greenMoments.resize(3, io::Matrix<Value>(jHybMatrix.size(), jHybMatrix.size()));
                    
                    for(std::size_t i = 0; i < jHybMatrix.size(); ++i) greenMoments[0](i, i)  = 1.;
                    
                    for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                        for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                        {
                            std::string entry = jHybMatrix(i)(j).string();
                            
                            if(!entry.empty()) {
                                jsx::value temp;
                                
                                linalg::mult<Value>('n', 'c', 1., jC(i), jParams("operators")(j), .0, temp);
                                linalg::mult<Value>('c', 'n', 1., jParams("operators")(j), jC(i), 1., temp);
                                greenMoments[1](i, j) += linalg::trace<Value>(jDensityMatrix, temp);
                                
                                linalg::mult<Value>('n', 'c', 1., jC(i), jC(j), .0, temp);
                                linalg::mult<Value>('c', 'n', 1., jC(j), jC(i), 1., temp);
                                greenMoments[2](i, j) += linalg::trace<Value>(jDensityMatrix, temp);
                            }
                        }
                    
                    
                    if(jParams.is("dyn"))
                    {
                        jsx::value jq = jParams("dyn")("quantum numbers");
                        jsx::value jDynFunctions = jParams("dyn")("functions");
                        jsx::value jDynMatrix = jParams("dyn")("matrix");
                        
                        std::vector<std::vector<double>> q(jq.size());
                        
                        std::vector<double> Qavg(jq.size(), .0);
                        std::vector<std::vector<double>> DQDQavg(jq.size(), std::vector<double>(jq.size(), .0));
                        std::vector<std::vector<double>> Dw0(jq.size(), std::vector<double>(jq.size(), .0));
                        std::vector<std::vector<double>> Dt0(jq.size(), std::vector<double>(jq.size(), .0));
                        
                        auto sectorProb = jsx::at<io::rvec>(jPartition("sector prob"));
                        jsx::value jDensityMatrixDyn = jsx::array_t(jq.size());
                        
                        if(jDynMatrix.size() != jq.size())
                            throw std::runtime_error("imp::Simple: matrix has wrong size !");
                        
                        for(int I = 0; I < jq.size(); ++I)
                        {
                            q[I] = jsx::at<io::rvec>(jq(I));

                            auto Q = ga::construct_sector_qn(jParams("hloc"), jq(I));
                            
                            for(int s = 0; s < jParams("hloc")("eigen values").size(); ++s) Qavg[I] += Q.at(s)*sectorProb.at(s);
                            
                            
                            if(jDynMatrix(I).size() != jq.size())
                                throw std::runtime_error("imp::Simple: matrix has wrong size !");
                            
                            for(int J = 0; J < jq.size(); ++J)
                            {
                                auto const& D = jsx::at<io::rvec>(jDynFunctions(jDynMatrix(I)(J).string()));
                                
                                Dw0[I][J] = D[0];
                                
                                for(std::size_t m = 0; m < D.size(); ++m) Dt0[I][J] += 2./beta*D[m];
                                
                                DQDQavg[I][J] = jsx::at<io::rvec>(jPartition("QQ")(I))[J];
                            }
                            
                            jDensityMatrixDyn(I) = partition::meas::read_density_matrix<Value>(jParams, jPartition("density matrix dyn")(I));
                        }

                        
                        for(int i = 0; i < jHybMatrix.size(); ++i)
                            for(int I = 0; I < jq.size(); ++I)
                                for(int K = 0; K < jq.size(); ++K)
                                    greenMoments[1](i, i) += -Qavg[K]*Dw0[K][I]*q[I][i];
                        
                        
                        for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                            for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                            {
                                std::string entry = jHybMatrix(i)(j).string();

                                if(!entry.empty())
                                {
                                    jsx::value temp;
                                    
                                    linalg::mult<Value>('n', 'c', 1., jParams("operators")(i), jC(j), .0, temp);
                                    linalg::mult<Value>('c', 'n', 1., jC(j), jParams("operators")(i), 1., temp);
                                    
                                    for(int I = 0; I < jq.size(); ++I)
                                        if(q[I][i] != .0) greenMoments[2](i, j) += -q[I][i]*linalg::trace<Value>(jDensityMatrixDyn(I), temp);
                                        
                                
                                    linalg::mult<Value>('n', 'c', 1., jC(i), jParams("operators")(j), .0, temp);
                                    linalg::mult<Value>('c', 'n', 1., jParams("operators")(j), jC(i), 1., temp);
                                    
                                    for(int J = 0; J < jq.size(); ++J)
                                        if(q[J][j] != .0) greenMoments[2](i, j) += -q[J][j]*linalg::trace<Value>(jDensityMatrixDyn(J), temp);
                                }
                            }

                        
                        for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                            for(int I = 0; I < jq.size(); ++I)
                                for(int J = 0; J < jq.size(); ++J)
                                    greenMoments[2](i, i) += q[I][i]*(-Dt0[I][J] + DQDQavg[I][J])*q[J][i];
                        
                    }
                    
                    greenMoments = func::get_function_matrix(func::get_function_entries(greenMoments, jHybMatrix), jHybMatrix);
                    
                    for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                        for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                            greenMoments[2](i, j) += hybMoments[0](i, j);
                    
                    return greenMoments;
                }
                
                template <typename Value>
                std::vector<io::Matrix<Value>> compute_self_moments(jsx::value const& jParams, std::vector<io::Matrix<Value>> const& hybMoments, std::vector<io::Matrix<Value>>& greenMoments) {
                    
                    double const mu = jParams("mu").real64();
                    auto const oneBody = jsx::at<io::Matrix<Value>>(jParams("hloc")("one body"));
                    jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
                    
                    std::vector<io::Matrix<Value>> selfMoments(2, io::Matrix<Value>(jHybMatrix.size(), jHybMatrix.size()));
                    
                    selfMoments.resize(2, io::Matrix<Value>(jHybMatrix.size(), jHybMatrix.size()));
                    
                    io::Matrix<Value> gm1gm1(jHybMatrix.size(), jHybMatrix.size());
                    linalg::mult<Value>('n', 'n', 1., greenMoments[1], greenMoments[1], .0, gm1gm1);
                    
                    for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                        for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                            selfMoments[0](i, j) += (i == j ? mu : .0) - oneBody(i, j) - greenMoments[1](i, j);
                            selfMoments[1](i, j) += greenMoments[2](i, j) - gm1gm1(i, j);
                        }
                    
                    selfMoments = get_function_matrix(get_function_entries(selfMoments, jHybMatrix), jHybMatrix);
                    
                    
                    return  selfMoments;
                    
                }
                
                template <typename Value>
                void add_green_tail(jsx::value const& jParams, iOmega const& iomega, io::Matrix<Value> const& oneBody, std::vector<io::cmat> const& hyb, std::vector<io::cmat> const& self, std::vector<io::cmat>& green){
                    
                    double const mu = jParams("mu").real64();
                    jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
                    
                    for(std::size_t n = green.size(); n < hyb.size(); ++n) {
                        io::cmat green_inv(jHybMatrix.size(), jHybMatrix.size());
                        
                        for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                            for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                                green_inv(i, j) = (i == j ? iomega(n) + mu : .0) - oneBody(i, j) - hyb[n](i, j) - self[n](i, j);
                        
                        green.push_back(linalg::inv(green_inv));
                    }
                }
                
                
                template<typename Value>
                void add_self_tail(jsx::value const& jHybMatrix, iOmega const& iomega, std::vector<io::cmat>& function, std::vector<io::Matrix<Value>> const& moments, std::size_t hybSize)
                {
                    io::Matrix<Value> alpha(jHybMatrix.size(), jHybMatrix.size());
                    io::Matrix<Value> gamma(jHybMatrix.size(), jHybMatrix.size());
                    
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
                                
                                alpha(i, j) = ut::to_val<Value>(-moments[1](i, j)/2.*(1./(selfPosAvg - moments[0](i, j)) + 1./(selfNegAvg - moments[0](i, j))));
                                gamma(i, j) = ut::to_val<Value>(iomegaAvg*moments[1](i, j)/2.*(-1./(selfPosAvg - moments[0](i, j)) + 1./(selfNegAvg - moments[0](i, j))) + iomegaAvg*iomegaAvg);
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
            }
        }
        
    }
    
}


#endif









