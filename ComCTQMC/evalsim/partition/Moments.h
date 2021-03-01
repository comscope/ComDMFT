#ifndef EVALSIM_PARTITION_MOMENTS_H
#define EVALSIM_PARTITION_MOMENTS_H


#include "ReadDensityMatrix.h"
#include "ReadHamiltonian.h"

#include "../../include/linalg/Operators.h"
#include "../../include/JsonX.h"
#include "../../include/io/Vector.h"
#include "../../include/io/Matrix.h"

namespace evalsim {
    
    namespace partition {
        
        
        template<typename Value>
        std::vector<io::Matrix<Value>> get_green_moments(jsx::value const& jParams, std::vector<io::Matrix<Value>> const& hybMoments, jsx::value const& jMeasurements, jsx::value const& jScalar)
        {
            double const beta = jParams("beta").real64();
            jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
            jsx::value jDensityMatrix = meas::read_density_matrix<Value>(jParams, jMeasurements("density matrix"));
            
            std::vector<io::Matrix<Value>> greenMoments;
            
            jsx::value jHamiltonian = get_hamiltonian<Value>(jParams);
            
            jsx::value jC = jsx::array_t(jHybMatrix.size());
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
            {
                linalg::mult<Value>('n', 'n',  1., jHamiltonian, jParams("operators")(i), .0, jC(i));
                linalg::mult<Value>('n', 'n', -1., jParams("operators")(i), jHamiltonian, 1., jC(i));
            }
            
            greenMoments.resize(3, io::Matrix<Value>(jHybMatrix.size(), jHybMatrix.size()));
            
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i) greenMoments[0](i, i)  = 1.;
            
            
            std::size_t const N = jHybMatrix.size();
            std::size_t const size = N*N;
            std::size_t const rank = mpi::rank();
            std::size_t const chunk = (size + mpi::number_of_workers() - 1)/mpi::number_of_workers();
            std::vector<Value> tmp1(size,0);
            std::vector<Value> tmp2(size,0);
            
            for(std::size_t index = rank*chunk; index < chunk*(rank + 1); ++index){
                if (index>=size) break;
                
                std::size_t const i = index%(N);
                std::size_t const j = index/(N);
                    
                std::string entry = jHybMatrix(i)(j).string();
                    
                if(!entry.empty()) {
                    jsx::value temp;
                        
                    linalg::mult<Value>('n', 'c', 1., jC(i), jParams("operators")(j), .0, temp);
                    linalg::mult<Value>('c', 'n', 1., jParams("operators")(j), jC(i), 1., temp);
                    tmp1[i*N+j] += linalg::trace<Value>(jDensityMatrix, temp);
                        
                    linalg::mult<Value>('n', 'c', 1., jC(i), jC(j), .0, temp);
                    linalg::mult<Value>('c', 'n', 1., jC(j), jC(i), 1., temp);
                    tmp2[i*N+j] += linalg::trace<Value>(jDensityMatrix, temp);
                }
            }
            
            mpi::reduce<mpi::op::sum>(tmp1, mpi::master);
            mpi::reduce<mpi::op::sum>(tmp2, mpi::master);
            
            for (std::size_t i=0; i<N; i++)
                for (std::size_t j=0; j<N; j++){
                    greenMoments[1](i,j) = tmp1[i*N+j];
                    greenMoments[2](i,j) = tmp2[i*N+j];
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
                
                auto sectorProb = jsx::at<io::rvec>(jMeasurements("sector prob"));
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
                        
                        DQDQavg[I][J] = jsx::at<io::rvec>(jMeasurements("QQ")(I))[J];
                    }
                    
                    jDensityMatrixDyn(I) = meas::read_density_matrix<Value>(jParams, jMeasurements("density matrix dyn")(I));
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
        
        
        template<typename Value>
        std::vector<io::Matrix<Value>> get_self_moments(jsx::value const& jParams, std::vector<io::Matrix<Value>> const& hybMoments, std::vector<io::Matrix<Value>> const& greenMoments)
        {
            double const mu = jParams("mu").real64();
            auto const& oneBody = jsx::at<io::Matrix<Value>>(jParams("hloc")("one body"));
            jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
            
            std::vector<io::Matrix<Value>> selfMoments(2, io::Matrix<Value>(jHybMatrix.size(), jHybMatrix.size()));
            
            io::Matrix<Value> gm1gm1(jHybMatrix.size(), jHybMatrix.size());
            linalg::mult<Value>('n', 'n', 1., greenMoments[1], greenMoments[1], .0, gm1gm1);
            
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                    selfMoments[0](i, j) += (i == j ? mu : .0) - oneBody(i, j) - greenMoments[1](i, j);
                    selfMoments[1](i, j) += greenMoments[2](i, j) - hybMoments[0](i, j) - gm1gm1(i, j);
                }
            
            selfMoments = func::get_function_matrix(func::get_function_entries(selfMoments, jHybMatrix), jHybMatrix);
            
            return selfMoments;
        }
        
    }
}


#endif









