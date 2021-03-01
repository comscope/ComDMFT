#ifndef EVALSIM_IMPRGREENTWO_MEASUREMENTS_H
#define EVALSIM_IMPRGREENTWO_MEASUREMENTS_H

#include <vector>
#include <complex>
#include <fstream>
#include <valarray>
#include <cstring>
#include <random>


#include "../../../../include/JsonX.h"
#include "../../../../include/io/Vector.h"
#include "../../../../include/io/Matrix.h"
#include "../../../../include/io/Tensor.h"

#include "Functions.h"
#include "Bessel.h"

namespace evalsim {
    
    namespace worm {
        
        typedef struct Fermion{
            
            static int const zero_freq = 1;
            
        } Fermion;
        typedef struct Boson{
            
            static int const zero_freq = 0;
            
        } Boson;
        
        
        
        namespace meas {
            
            namespace impl {
                
                io::cvec read_matsubara_function(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jFunction, bool scale){
                    double const beta = jParams("beta").real64();
                    
                    io::cvec function = jsx::at<io::cvec>(jFunction);
                    
                    if (scale) for (auto& x : function) x*=-1./beta;
                    
                    return function;
                }
                
                io::cvec collapse_antisymmetric_function(io::cvec const& function){
                    
                    io::cvec as_function(function.size()/2);
                    
                    std::size_t start = function.size()/2;
                    
                    auto it_forward = function.begin()+start;
                    auto it_backward = function.begin()+start-1;
                    
                    for (auto& x : as_function){
                        x = 0.5 * ( *it_forward + std::conj(*it_backward) );
                        --it_backward; ++it_forward;
                    }
                    
                    return as_function;
                }
                
                
                
                template <typename T, typename ...> struct Read_legendre_function;
                
                template<typename T>
                struct Read_legendre_function<T>{
                  
                    Read_legendre_function(jsx::value const& jParams, jsx::value const& jPartition){ throw std::runtime_error("Error: Read_legendre_function requires specialization\n");}
                    
                    io::cvec operator()(jsx::value const& jFunction){
                        
                        return io::cvec();
                    }
                };
                
                template<typename T>
                struct Read_legendre_function<T,Fermion,Fermion,Boson>{
                    
                    Read_legendre_function(jsx::value const& jParams, jsx::value const& jPartition) :
                    nMatGB_((std::is_same<T,double>::value ? 1 : 2) * jPartition("boson cutoff").int64()),
                    nMatGF_(2*(jPartition.is("matsubara cutoff") ? jPartition("matsubara cutoff").int64() : 50 )),
                    nPol_(jPartition("fermion cutoff").int64()),
                    beta_(jParams("beta").real64()),
                    Tnl_(nPol_,nMatGF_){}
                                       
                    io::cvec operator()(jsx::value const& jFunction){
                        
                        io::cvec legendre_function = jsx::at<io::cvec>(jFunction);
                        
                        io::cvec function(nMatGF_*nMatGF_*nMatGB_,0);
                        
                        for (int l=0; l<nPol_; l++)
                            for (int m=0; m<nPol_; m++){
                                
                                int const s = m%2 ? 1 : -1;
                                double const adj=-s*std::sqrt(2*l+1.)*std::sqrt(2*m+1.)/beta_;
                                
                                for (int k=0; k<nMatGB_; k++){
                                    int n = l + m*nPol_ + k*nPol_*nPol_;
                                    legendre_function[n]*=adj;
                                    
                                    for (int i=0; i<nMatGF_; i++){
                                        int const i_=i-nMatGF_/2;
                                        for (int j=0; j<nMatGF_; j++){
                                            int const j_=j-nMatGF_/2;
                                            
                                            function[i + j*nMatGF_ + k*nMatGF_*nMatGF_] += legendre_function[n]*Tnl_(i_,l)*std::conj(Tnl_(j_,m));
                                            
                                        }
                                    }
                                }
                            }
                        
                        
                        return function;
                        
                    }
                    
                private:
                    int const nMatGB_,nMatGF_,nPol_;
                    double const beta_;
                    
                    be::Transform<T,Fermion> Tnl_;
                    
                };
                
                
                template<typename T>
                struct Read_legendre_function<T,Fermion,Boson>{
                    
                    Read_legendre_function(jsx::value const& jParams, jsx::value const& jPartition) :
                    nMatGB_((std::is_same<T,double>::value ? 1 : 2) * jPartition("boson cutoff").int64()),
                    nMatGF_((std::is_same<T,double>::value ? 1 : 2) * (jPartition.is("matsubara cutoff") ? jPartition("matsubara cutoff").int64() : 50 ) ),
                    nPol_(jPartition("fermion cutoff").int64()),
                    beta_(jParams("beta").real64()),
                    Tnl_(nPol_,nMatGF_){}
                    
                    io::cvec operator()(jsx::value const& jFunction){
                        
                        io::cvec legendre_function = jsx::at<io::cvec>(jFunction);
                        
                        io::cvec function(nMatGF_*nMatGB_,0);
                        
                        for (int m=0; m<nMatGB_; m++)
                            for (int l=0; l<nPol_; l++){
                            
                                legendre_function[m*nPol_+l]*=-std::sqrt(2*l+1.)/beta_;
                                
                                for (int n=0; n<nMatGF_; n++){
                                    int n_ = n - (std::is_same<T,double>::value ? 0 : nMatGF_/2);
                            
                                    function[m*nMatGF_+n] += legendre_function[m*nPol_+l]*Tnl_(n_,l);
                                }
                            
                        }
                        
                        return function;
                        
                    }
                    
                private:
                    int const nMatGB_,nMatGF_,nPol_;
                    double const beta_;
                    
                    be::Transform<T,Fermion> Tnl_;
                    
                };
                
                
                template<typename T, typename Time>
                struct Read_legendre_function<T,Time>{
                    
                    Read_legendre_function(jsx::value const& jParams, jsx::value const& jPartition) :
                    nMatGF_((std::is_same<T,double>::value ? 1 : 2) * (jPartition.is("matsubara cutoff") ? jPartition("matsubara cutoff").int64() : 50 )),
                    nPol_(jPartition("cutoff").int64()),
                    beta_(jParams("beta").real64()),
                    Tnl_(nPol_,nMatGF_){
                        
                        if (std::is_same<Time,Boson>::value) {throw std::runtime_error("Legendre transform is not implemented for bosonic times\n");}
                        //if (std::is_same<Time,Boson>::value)
                        //    for (int i=0; i<nMatGF_; i++) std::cout << i << " " << Tnl_(i,0) << "\n";
                    }
                    
                    io::cvec operator()(jsx::value const& jFunction){
                        
                        io::Vector<T> legendre_function = jsx::at<io::Vector<T>>(jFunction);
                        
                        io::cvec function(nMatGF_,0);
                        auto const shift = (std::is_same<T,double>::value ? 0 : nMatGF_/2);
                        
                        for (int l=0; l<nPol_; l++){
                            
                            legendre_function[l]*=-std::sqrt(2.0*l+1)/beta_;
                            
                            for (int n=0; n<nMatGF_; n++){
                                int const n_ = n - shift;
                                function[n] += legendre_function[l]*Tnl_(n_,l);
                            }
                            
                        }
                        
                        return function;
                    }
                    
                private:
                    int const nMatGF_,nPol_;
                    double const beta_;
                    
                    be::Transform<T,Time> Tnl_;
                    
                };
                
                
                template<typename Value, typename ... Args>
                inline io::cvec read_function(jsx::value const& jFunction, jsx::value const& jParams, jsx::value const& jPartition, int hybSize, bool scale, bool force_matsubara = false) {
                    if(force_matsubara or (jPartition.is("basis") ? jPartition("basis").string() == "matsubara" : true)) {
                        return read_matsubara_function(jParams,jPartition, jFunction, scale);
                    }else if(jPartition.is("basis") ? jPartition("basis").string() == "legendre" : true) {
                         Read_legendre_function<double,Args ...> rlf(jParams,jPartition);
                        return rlf(jFunction);
                    }else
                        throw std::runtime_error("read_function: unknown basis option");
                    
                    throw std::runtime_error("read_function: invalid format");
                };
                
                
            }
            
            
            template<typename Value, typename ... Args>
            std::vector<io::ctens> read_tensor_functions(jsx::value const& jFunctions, jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jHybMatrix, int hybSize)
            {
                std::map<std::string, io::cvec> functions;
                
                for(auto const& function : jFunctions.object()) {
                    auto temp = impl::read_function<Value, Args ... >(function.second, jParams, jPartition, hybSize, true);
                    
                    functions[function.first] = temp;
                }
                
                std::vector<io::ctens> tensors = func::get_function_tensor(functions, jPartition, jHybMatrix);
                
                return tensors;
            };
            
            
            template<typename Value,typename ... Args>
            std::vector<io::cmat> read_matrix_functions(jsx::value const& jFunctions, jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jHybMatrix, int hybSize)
            {
                std::map<std::string, io::cvec> functions;
                
                for(auto const& function : jFunctions.object()) {
                    auto temp = impl::read_function<Value,Args ... >(function.second, jParams, jPartition, hybSize, true, false);
                    if (std::is_same<Value,ut::complex>::value) temp = impl::collapse_antisymmetric_function(temp);
                    
                    functions[function.first] = temp;
                }
                
                std::vector<io::cmat> matrix = func::get_function_matrix(functions, jHybMatrix);
                
                return matrix;
            };

            
            
            template<typename Value>
            std::vector<io::cmat>  read_matrix_functions_from_obs(jsx::value const& jFunctions, jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jHybMatrix, int hybSize){
                
                std::map<std::string, io::cvec> functions;
                
                for(auto const& function : jFunctions.object()) {
                    auto temp = impl::read_function<Value>(function.second("function"), jParams, jPartition, hybSize, false, true);
                    
                    functions[function.first] = temp;
                }
                
                std::vector<io::cmat> matrix = func::get_function_matrix(functions, jHybMatrix);
                
                return matrix;
                
            }
            
            template<typename Value>
            std::vector<io::ctens>  read_tensor_functions_from_obs(jsx::value const& jFunctions, jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jHybMatrix, int hybSize){
                
                std::map<std::string, io::cvec> functions;
                
                for(auto const& function : jFunctions.object()) {
                    auto temp = impl::read_function<Value>(function.second("function"), jParams, jPartition, hybSize, false, true);
                    
                    functions[function.first] = temp;
                }
                
                std::vector<io::ctens> tensor = func::get_function_tensor(functions, jPartition, jHybMatrix);
                
                return tensor;
                
            }
            
            
        }
        
    }
    
}

#endif
