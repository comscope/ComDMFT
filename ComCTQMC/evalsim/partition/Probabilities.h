#ifndef EVALSIM_PARTITION_PROBABILITIES_H
#define EVALSIM_PARTITION_PROBABILITIES_H


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
        
        
        
        
        template<typename Value>
        jsx::value get_probabilities(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements)
        {
            jsx::value jProbabilities;  
            
            mpi::cout << "Reading impurity probabilities ... " << std::flush;
            
            jsx::value jHamiltonianEff = get_effective_hamiltonian<Value>(jParams);
            jsx::value jDensityMatrix = meas::read_density_matrix<Value>(jParams, jMeasurements("density matrix"));
            
            
            //Collect list of probabilities that have survived the criterion
            std::vector<std::string> surviving;
            for(int qn = 0; qn < jPartition("probabilities").size(); ++qn) {
                std::string const name = jPartition("probabilities")(qn).string();
                if(name == "energy" || jPartition("quantum numbers").is(name) || jPartition("observables").is(name))
                    surviving.push_back(name);
            }
            
            jsx::array_t temp;
            for(auto& entry : surviving) {
                temp.push_back(jsx::string_t(entry));
            }
            jProbabilities["surviving"] = std::move(temp);
            
            
            
            std::vector<std::vector<double>> data;
            for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector)
                for(int i = 0; i < jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size(); ++i) {
                    std::vector<double> temp(surviving.size() + 1);
                    temp.back() = std::abs(jsx::at<io::Matrix<Value>>(jDensityMatrix(sector)("matrix"))(i, i));     // take abs(real) value ?
                    data.push_back(temp);
                }
            
            for(int qn = 0; qn < surviving.size(); ++qn) {
                std::string const name = surviving[qn];
                
                if(jPartition("quantum numbers").is(name)) {
                    auto Qn = ga::construct_sector_qn(jParams("hloc"), jPartition("quantum numbers")(name));
                    
                    int index = 0;
                    for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector)
                        for(int i = 0; i < jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size(); ++i)
                            data[index++][qn] = truncate(Qn.at(sector), 8);
                    
                } else if(jPartition("observables").is(name)) {
                    jsx::value const& jObservable = jPartition("observables")(name);
                    
                    int index = 0;
                    for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector)
                        for(int i = 0; i < jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size(); ++i)
                            data[index++][qn]= truncate(ut::real(jsx::at<io::Matrix<Value>>(jObservable(sector)("matrix"))(i, i)), 8);
                    
                } else if(name == "energy") {
                    int index = 0;
                    for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector)
                        for(int i = 0; i < jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size(); ++i)
                            data[index++][qn] = truncate(ut::real(jsx::at<io::Matrix<Value>>(jHamiltonianEff(sector)("matrix"))(i, i)), 8);
                    
                }
            }
            
            
            if(surviving.size()) std::sort(data.begin(), data.end(), VecLess(surviving.size() + 1, 0, surviving.size()));
            
            
            temp.clear();
            for(auto& entry : data) {
                temp.push_back(jsx::array_t(entry.begin(), entry.end()));
                temp.back().array().back() = io::rvec{{entry.back()}};
            }
            jProbabilities["quantum numbers"] = std::move(temp);
            
            
            
            {
                jsx::value jTransformation = jParams("hloc")("transformation"), jTemp;
                
                linalg::mult<Value>('n', 'c', 1., jDensityMatrix, jTransformation, .0, jTemp);
                linalg::mult<Value>('n', 'n', 1., jTransformation, jTemp, .0, jDensityMatrix);
            }
            
            jsx::value jOccupationStates = ga::construct_occupation_states<Value>(jParams("hloc"));
            
            data.clear(); int index = 0;
            for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector)
                for(int i = 0; i < jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size(); ++i) {
                    io::rvec temp = jsx::at<io::rvec>(jOccupationStates(index++));
                    temp.push_back(std::abs(jsx::at<io::Matrix<Value>>(jDensityMatrix(sector)("matrix"))(i, i)));       // take abs(real) value ?
                    data.push_back(temp);
                }
            
            jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
            
            std::sort(data.begin(), data.end(), VecLess(jHybMatrix.size() + 1, jHybMatrix.size(), jHybMatrix.size() + 1));
            
            temp.clear();
            for(auto& entry : data) {
                temp.push_back(jsx::array_t(entry.begin(), entry.end()));
                temp.back().array().back() = io::rvec{{entry.back()}};
            }
            jProbabilities["occupation numbers"] = std::move(temp);
            
            mpi::cout << "Ok" << std::endl;
            
            return jProbabilities;
        }
   
    }
    
}


#endif









