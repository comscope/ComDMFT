#ifndef EVALSIM_INCLUDE_FUNCTIONS
#define EVALSIM_INCLUDE_FUNCTIONS

#include <tuple>
#include <vector>
#include <string>

#include "../../../../include/measurements/Measurements.h"
#include "../../../../include/io/Vector.h"
#include "../../../../include/io/Matrix.h"
#include "../../../../include/io/Tensor.h"
#include "../../../../include/atomic/Generate.h"
#include "../../../../ctqmc/include/bath/Hyb.h"
#include "../../../../ctqmc/include/Utilities.h"


namespace evalsim {
    
    namespace worm {
        
        namespace func {
            
            
            //Works for either partition measurements (labeled by their hybridisation)
            //or worm measurements (labeld by the operators c_i c_j^dag)
            template<typename Value>
            io::Matrix<Value> get_matrix(std::map<std::string, Value> const& entries, jsx::value const& jMatrix)
            {
                io::Matrix<Value> matrix(jMatrix.size(), jMatrix.size());
                
                for(std::size_t i = 0; i < jMatrix.size(); ++i)
                    for(std::size_t j = 0; j < jMatrix.size(); ++j) {
                        auto const entry = std::to_string(2*i)+"_"+std::to_string(2*j+1);
                        auto const entry_partition = jMatrix(i)(j).string();
                        
                        auto it = entries.find(entry);
                        if(it != entries.end()) matrix(i, j) = entries.at(entry);
                        
                        it = entries.find(entry_partition);
                        if(it != entries.end()) {matrix(i, j) = entries.at(entry_partition);}
                        
                    }
                
                //make symmetric if real
                
                return matrix;
            }
            
            
            template<typename Value>
            inline std::vector<io::Matrix<Value>> get_function_matrix(std::map<std::string, io::Vector<Value>> const& functions, jsx::value const& jMatrix)
            {
                std::size_t size = std::numeric_limits<std::size_t>::max();
                
                for(auto const& function : functions)
                    size = std::min(size, function.second.size());
                
                std::vector<io::Matrix<Value>> functionMatrix;
                for(std::size_t n = 0; n < size; ++n) {
                    std::map<std::string, Value> entries;
                    
                    for(auto& function : functions)
                        entries[function.first] = function.second[n];
                    
                    
                    functionMatrix.push_back(get_matrix(entries, jMatrix));
                }
                
                return functionMatrix;
            }
           
            template<typename Value>
            io::Tensor<Value> get_tensor(std::map<std::string, Value> const& entries, jsx::value const& jWorm, jsx::value const& jMatrix)
            {
                io::Tensor<Value> tensor(jMatrix.size(), jMatrix.size(), jMatrix.size(), jMatrix.size());
                
                for(std::size_t i = 0; i < jMatrix.size(); ++i)
                    for(std::size_t j = 0; j < jMatrix.size(); ++j)
                        for(std::size_t k = 0; k < jMatrix.size(); ++k)
                            for(std::size_t l = 0; l < jMatrix.size(); ++l) {
                                
                                //Check for any possible compbination of creation or annihilation operators
                                //Store the resulting entry in the tensor (would be nice to avoid...)
                                //So that on output we can know which combination was measured [see get_entries(io::Tensor<Value> ... )]
                                //This could alternately be specified by worm, but this adds quite a bit of code
                                for (int i_dagg = 0; i_dagg < 2; i_dagg++)
                                for (int j_dagg = 0; j_dagg < 2; j_dagg++)
                                for (int k_dagg = 0; k_dagg < 2; k_dagg++)
                                for (int l_dagg = 0; l_dagg < 2; l_dagg++){
                                    auto const entry = std::to_string(2*i+i_dagg)+"_"+std::to_string(2*j+j_dagg)+"_"+std::to_string(2*k+k_dagg)+"_"+std::to_string(2*l+l_dagg);
                                    if(entries.find(entry) != entries.end()){ tensor.emplace(i, j, k, l, entry, entries.at(entry));}
                                }
                                
                            }
                
                return tensor;
            }
            
            
            template<typename Value>
            inline std::vector<io::Tensor<Value>> get_function_tensor(std::map<std::string, io::Vector<Value>> const& functions, jsx::value const& jParams, jsx::value const& jMatrix)
            {
                std::size_t size = std::numeric_limits<std::size_t>::max();
                
                for(auto const& function : functions)
                    size = std::min(size, function.second.size());
                
                std::vector<io::Tensor<Value>> functionTensor;
                for(std::size_t n = 0; n < size; ++n) {
                    std::map<std::string, Value> entries;
                    
                    for(auto& function : functions){
                        entries[function.first] = function.second[n];
                    }
                    
                    functionTensor.push_back(get_tensor<Value>(entries, jParams, jMatrix));
                }
                
                return functionTensor;
            }
            
            
            template<typename Value>
            inline std::map<std::string, Value> get_entries(io::Matrix<Value> const& matrix, jsx::value const& jMatrix)
            {
                std::map<std::string, std::pair<int, Value>> temp;
                
                for(std::size_t i = 0; i < jMatrix.size(); ++i)
                    for(std::size_t j = 0; j < jMatrix.size(); ++j) {
                        //auto const entry = std::to_string(2*i)+"_"+std::to_string(2*j+1);
                        auto const entry = jMatrix(i)(j).string();
                        
                        if (std::abs(matrix(i,j))) {
                            if(!temp.count(entry))
                                temp[entry] = std::pair<int, Value>(0, .0);
                            
                            temp.at(entry).first  += 1;
                            temp.at(entry).second += matrix(i, j);
                        }
                    }
                
                std::map<std::string, Value> entries;
                for(auto const& entry : temp)
                    entries[entry.first] = entry.second.second/Value(entry.second.first);
                
                return entries;
            }
            
            
            template<typename Value>
            inline void add_tensor_result(Value const value, std::vector<int> const rearranged_indices, std::map<std::string, std::pair<int, Value>> & entries_map){
            
                std::string entry = std::to_string(rearranged_indices[0])+"_"+std::to_string(rearranged_indices[1])+"_"+std::to_string(rearranged_indices[2])+"_"+std::to_string(rearranged_indices[3]);
                
                if(!entries_map.count(entry))
                entries_map[entry] = std::pair<int, Value>(0, .0);
                
                entries_map.at(entry).first  += 1;
                entries_map.at(entry).second += value;
                
            }
            
            
            
            template<typename Value>
            inline std::map<std::string, Value> get_entries(io::Tensor<Value> const& tensor, jsx::value const& jMatrix)
            {
                std::map<std::string, std::pair<int, Value>> temp;
                
                //get entry name from tensor
                //It would be nice to avoid complicating the tensor class, but there is not a great way to know which entry
                for(std::size_t i = 0; i < jMatrix.size(); ++i)
                    for(std::size_t j = 0; j < jMatrix.size(); ++j)
                        for(std::size_t k = 0; k < jMatrix.size(); ++k)
                            for(std::size_t l = 0; l < jMatrix.size(); ++l) {
                                auto const entry = (tensor.entry(i, j, k, l));
                                
                                if (entry!=""){
                                    
                                    //Split entry into indices
                                    if(!temp.count(entry))
                                        temp[entry] = std::pair<int, Value>(0, .0);
                                    
                                    temp.at(entry).first  += 1;
                                    temp.at(entry).second += tensor(i,j,k,l);
                                    
                                }
                                
                            }
                
                
                std::map<std::string, Value> entries;
                for(auto const& entry : temp)
                    entries[entry.first] = entry.second.second/Value(entry.second.first);
                
                return entries;
            }
            
            
            template<typename Value>
            inline std::map<std::string, io::Vector<Value>> get_function_entries(std::vector<io::Matrix<Value>> const& functionMatrix, jsx::value const& jMatrix)
            {
                std::map<std::string, io::Vector<Value>> functionEntries;
                
                for(auto const& matrix : functionMatrix) {
                    std::map<std::string, Value> entries = get_entries<Value>(matrix, jMatrix);
                    for(auto const& entry : entries){
                        functionEntries[entry.first].push_back(entry.second);
                    }
                }
                
                return functionEntries;
            }
            
            template<typename Value>
            inline std::map<std::string, io::Vector<Value>> get_function_entries(std::vector<io::Tensor<Value>> const& functionTensor, jsx::value const& jMatrix)
            {
                std::map<std::string, io::Vector<Value>> functionEntries;
                
                for(auto const& matrix : functionTensor) {
                    std::map<std::string, Value> entries = get_entries(matrix, jMatrix);
                    for(auto const& entry : entries)
                        functionEntries[entry.first].push_back(entry.second);
                }
                
                return functionEntries;
            }
            
            
            template<typename Value>
            jsx::value write_functions(jsx::value const& jParams, jsx::value const& jHybMatrix, std::vector<io::cmat> const& functionsMatrix, std::vector<io::Matrix<Value>> const& momentsMatrix) {
                jsx::value jFunction;
                
                std::map<std::string, io::cvec> functions = get_function_entries(functionsMatrix, jHybMatrix);
                for(auto& function : functions){
                    jFunction[function.first]["function"] = std::move(function.second);
                }
                
                std::map<std::string, io::Vector<Value>> moments = get_function_entries(momentsMatrix, jHybMatrix);
                
                for(auto& moment : moments) jFunction[moment.first]["moments"] = moment.second;
                
                return jFunction;
            }
            
            template<typename Value>
            jsx::value write_functions(jsx::value const& jParams, jsx::value const& jHybMatrix, std::vector<io::cmat> const& functionsMatrix) {
                jsx::value jFunction;
                
                std::map<std::string, io::cvec> functions = get_function_entries(functionsMatrix, jHybMatrix);
                
                for(auto& function : functions){
                    jFunction[function.first]["function"] = std::move(function.second);
                }
                return jFunction;
            }
            

            template<typename Value>
            jsx::value write_functions(jsx::value const& jParams, jsx::value const& jHybMatrix, std::vector<io::ctens> const& functionsMatrix) {
                jsx::value jFunction;
                
                std::map<std::string, io::cvec> functions = get_function_entries(functionsMatrix, jHybMatrix);
                for(auto& function : functions){
                    double t=0; for (auto x : function.second) t+=std::abs(x);
                    if (t and !std::isnan(t)) jFunction[function.first]["function"] = std::move(function.second);
                }
                
                return jFunction;
            }
            
            
            template<typename Value>
            std::vector<io::Matrix<Value>> get_hybridisation_moment(std::map<std::string, io::cvec> const& functions, jsx::value const& jParams, jsx::value const& jHybMatrix)
            {
                std::map<std::string, Value> moments;
                
                for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                    for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                        if(jHybMatrix(i)(j).string() != "") {
                            if(jHybMatrix(j)(i).string() == "")
                                throw std::runtime_error("Hyb: invalid hybridisation matrix.");
                            
                            std::string const entry = jHybMatrix(i)(j).string();
                            std::string const entryTransp = jHybMatrix(j)(i).string();
                            
                            if(!moments.count(entry))
                                moments[entry] = bath::Fit<Value>(jParams("beta").real64(), functions.at(entry), functions.at(entryTransp)).moment();
                        }
                
                return  { get_matrix(moments, jHybMatrix) };
            }
            
            template<typename Value>
            std::tuple<std::vector<io::cmat>, std::vector<io::Matrix<Value>>> get_hybridisation(jsx::value const& jParams, jsx::value const& jHybMatrix)
            {
                
                jsx::value jFunctions = mpi::read(jParams("hybridisation")("functions").string());
                
                std::map<std::string, io::cvec> functions;
                for(auto& function : jFunctions.object())
                    functions[function.first] = jsx::at<io::cvec>(function.second);
                
                return std::make_tuple(get_function_matrix(functions, jHybMatrix),
                                       get_hybridisation_moment<Value>(functions, jParams, jHybMatrix));  // checks if hybridisation is compatible with Value
            }
            
            
        }
        
    }
    
}

#endif //EVALSIM










