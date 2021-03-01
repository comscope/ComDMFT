#ifndef EVALSIM_WORM_FULL_VERTEX_KERNELS
#define EVALSIM_WORM_FULL_VERTEX_KERNELS

#include "../../include/JsonX.h"
#include "../../ctqmc/include/config/Worms.h"
#include "../../ctqmc/include/Params.h"
#include "../../ctqmc/include/impurity/Tensor.h"

#include "include/functions/Functions.h"
#include "include/functions/Measurements.h"
#include "include/functions/Utilities.h"

namespace evalsim {
    
    namespace worm {
        
        struct Kernels{
           static std::string const name;
        };
        std::string const Kernels::name = "kernels";
        
        template<typename Value>
        struct InteractionTensor {
            InteractionTensor() = delete;
            InteractionTensor(imp::Tensor<Value> const& tensor, int N) :
            N_(N),
            tensor_(N*N*N*N) {
                
                for (int i=0; i<N; i++)
                for (int j=0; j<N; j++)
                for (int k=0; k<N; k++)
                for (int l=0; l<N; l++)
                    tensor_[N_*N_*N_*i + N_*N_*j + N_*k + l] = tensor(i,j,l,k);
                
            };
            
            Value operator()(int i, int j, int k, int l) const {
                return tensor_[N_*N_*N_*i + N_*N_*j + N_*k + l];
            };
            
        private:
            int const N_;
            io::Vector<Value> tensor_;
            
        };
        
        
        std::vector<io::ctens> rearrange_susc_ph(std::vector<io::ctens> const& susc_ph){
            std::vector<io::ctens> r(susc_ph.size(),io::ctens(susc_ph[0].I(),susc_ph[0].I(),susc_ph[0].I(),susc_ph[0].I()));
            for (auto const& ijkl : susc_ph[0].ijkl()){

                auto const i=ijkl[1];
                auto const j=ijkl[2];
                auto const k=ijkl[3];
                auto const l=ijkl[4];
                
                //j and i are swapped
                //k and l are swapped
                std::string const entry = std::to_string(2*j)+"_"+std::to_string(2*i+1)+"_"+std::to_string(2*l)+"_"+std::to_string(2*k+1);
                
                for (int n=0; n<susc_ph.size(); n++)
                    r[n].emplace(j,i,l,k, entry, -susc_ph[n].at(i,j,k,l));
                
            }

            return r;
        }
    
        std::vector<io::ctens> rearrange_susc_tph(std::vector<io::ctens> const& susc_ph){
            std::vector<io::ctens> r(susc_ph.size(),io::ctens(susc_ph[0].I(),susc_ph[0].I(),susc_ph[0].I(),susc_ph[0].I()));
            
            for (auto const& ijkl : susc_ph[0].ijkl()){

                auto const i=ijkl[1];
                auto const j=ijkl[2];
                auto const k=ijkl[3];
                auto const l=ijkl[4];
                
                //j and l are swapped
                std::string const entry = std::to_string(2*i)+"_"+std::to_string(2*l+1)+"_"+std::to_string(2*k)+"_"+std::to_string(2*j+1);
                
                for (int n=0; n<susc_ph.size(); n++)
                    r[n].emplace(i,l,k,j, entry, -susc_ph[n].at(i,j,k,l));
                
            }

            return r;
            
        }
    
        std::vector<io::ctens> rearrange_susc_pp(std::vector<io::ctens> const& susc_pp){
            std::vector<io::ctens> r(susc_pp.size(),io::ctens(susc_pp[0].I(),susc_pp[0].I(),susc_pp[0].I(),susc_pp[0].I()));
        
            for (auto const& ijkl : susc_pp[0].ijkl()){

                auto const i=ijkl[1];
                auto const j=ijkl[2];
                auto const k=ijkl[3];
                auto const l=ijkl[4];
                
                //j and k are swapped
                std::string const entry = std::to_string(2*i)+"_"+std::to_string(2*k+1)+"_"+std::to_string(2*j)+"_"+std::to_string(2*l+1);
                
                for (int n=0; n<susc_pp.size(); n++)
                    r[n].emplace(i,k,j,l, entry, -susc_pp[n].at(i,j,k,l));
                
            }

            return r;

        }
    
    
        std::vector<io::ctens> rearrange_hedin_ph(std::vector<io::ctens> const& hedin_ph){
            std::vector<io::ctens> r(hedin_ph.size(),io::ctens(hedin_ph[0].I(),hedin_ph[0].I(),hedin_ph[0].I(),hedin_ph[0].I()));
            
            for (auto const& ijkl : hedin_ph[0].ijkl()){

                auto const i=ijkl[1];
                auto const j=ijkl[2];
                auto const k=ijkl[3];
                auto const l=ijkl[4];
                
                //k and l are swapped
                std::string const entry = std::to_string(2*i)+"_"+std::to_string(2*j+1)+"_"+std::to_string(2*l)+"_"+std::to_string(2*k+1);

                for (int n=0; n<hedin_ph.size(); n++){
                    r[n].emplace(i,j,l,k, entry, hedin_ph[n].at(i,j,k,l));
                }
                
            }

            return r;

        }
    
        std::vector<io::ctens> rearrange_hedin_tph(std::vector<io::ctens> const& hedin_ph){
            std::vector<io::ctens> r(hedin_ph.size(),io::ctens(hedin_ph[0].I(),hedin_ph[0].I(),hedin_ph[0].I(),hedin_ph[0].I()));
        
            for (auto const& ijkl : hedin_ph[0].ijkl()){

                auto const i=ijkl[1];
                auto const j=ijkl[2];
                auto const k=ijkl[3];
                auto const l=ijkl[4];
                
                //l and j are swapped
                std::string const entry = std::to_string(2*i)+"_"+std::to_string(2*l+1)+"_"+std::to_string(2*k)+"_"+std::to_string(2*j+1);
                
                for (int n=0; n<hedin_ph.size(); n++)
                    r[n].emplace(i,l,k,j, entry, hedin_ph[n].at(i,j,k,l));
                
            }

            return r;
        }
    
        std::vector<io::ctens> rearrange_hedin_pp(std::vector<io::ctens> const& hedin_pp){
            std::vector<io::ctens> r(hedin_pp.size(),io::ctens(hedin_pp[0].I(),hedin_pp[0].I(),hedin_pp[0].I(),hedin_pp[0].I()));

            for (auto const& ijkl : hedin_pp[0].ijkl()){

                auto const i=ijkl[1];
                auto const j=ijkl[2];
                auto const k=ijkl[3];
                auto const l=ijkl[4];
                
                //k and j are swapped
                std::string const entry = std::to_string(2*i)+"_"+std::to_string(2*k+1)+"_"+std::to_string(2*j)+"_"+std::to_string(2*l+1);
                                   
                for (int n=0; n<hedin_pp.size(); n++)
                    r[n].emplace(i,k,j,l, entry, -hedin_pp[n].at(i,j,k,l));
                
                
            }

            return r;
        }
    
        template<typename Value>
        jsx::value evaluateKernels(jsx::value jParams, jsx::value const& jObservables) {
            
            
            std::cout << "Evaluating asymptotic vertex kernels" << std::endl;
            
            if (!jObservables.is(cfg::susc_ph::Worm::name()) or
                !jObservables.is(cfg::susc_pp::Worm::name()) or
                !jObservables.is(cfg::hedin_ph_imprsum::Worm::name()) or
                !jObservables.is(cfg::hedin_pp_imprsum::Worm::name())){
                
                std::cout << jObservables.is(cfg::susc_ph::Worm::name()) << std::endl;
                std::cout << jObservables.is(cfg::susc_pp::Worm::name()) << std::endl;
                std::cout << jObservables.is(cfg::hedin_ph_imprsum::Worm::name()) << std::endl;
                std::cout << jObservables.is(cfg::hedin_pp_imprsum::Worm::name()) << std::endl;
                
                throw std::runtime_error("Vertex kernel evaluation require ctqmc measurement of susc_ph/pp and hedin_ph/pp (imprsum versions)");
                
                return jsx::null_t();
            }
            
            jsx::value jWorm = jParams(cfg::hedin_ph_imprsum::Worm::name());
            
            //double const beta = jParams("beta").real64();
            int const pos_and_neg_freq = (std::is_same<Value,double>::value ? 1 : 2);
            int const nMatGB = pos_and_neg_freq*jWorm("boson cutoff").int64()-pos_and_neg_freq+1;
            int const nMatGF = 2*(jWorm("basis").string() == "matsubara" ? jWorm("fermion cutoff").int64() : ( jWorm.is("fermion cutoff") ? jWorm("fermion cutoff").int64() : 50 ));
            
            func::OmegaMap omega_f(nMatGF,false,true);
            func::OmegaMap omega_b(nMatGB,true,!std::is_same<Value,double>::value);
            
            int const nMatGB_kernel = 2*jWorm("boson cutoff").int64()-1;
            int const nMatGF_kernel = 2*(jWorm("basis").string() == "matsubara" ? jWorm("fermion cutoff").int64() : ( jWorm.is("fermion cutoff") ? jWorm("fermion cutoff").int64() : 50 ));
            
            func::OmegaMap omega_f_kernel(nMatGF_kernel,false,true);
            func::OmegaMap omega_b_kernel(nMatGB_kernel,true,true);
            
            
            jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
            std::vector<io::cmat> hyb; std::vector<io::Matrix<Value>> hybMoments;
            std::tie(hyb, hybMoments) = partition::func::get_hybridisation<Value>(jParams);
            
            //Read in green functions and susc and hedin susceptibilities
            std::cout << "Reading in Hedin and Susc Green's functions ... " << std::flush;
            
            std::vector<io::cmat> green_pos = meas::read_matrix_functions_from_obs<ut::complex>(jObservables(cfg::partition::Worm::name())("green"), jParams, jParams(cfg::partition::Worm::name()), jHybMatrix, hyb.size());
            auto const green = func::green_function_on_full_axis(green_pos);
            func::OmegaMap greenOM(green.size(),false,true);
            
            
            auto const susc_ph = rearrange_susc_ph( meas::read_tensor_functions_from_obs<ut::complex>(jObservables(cfg::susc_ph::Worm::name())("susceptibility"), jParams, jParams(cfg::susc_ph::Worm::name()), jHybMatrix, hyb.size()));
            auto const susc_pp = rearrange_susc_pp( meas::read_tensor_functions_from_obs<ut::complex>(jObservables(cfg::susc_pp::Worm::name())("susceptibility"), jParams, jParams(cfg::susc_ph::Worm::name()), jHybMatrix, hyb.size()));
            auto const susc_tph = rearrange_susc_tph(susc_ph);
            
            auto const hedin_ph = rearrange_hedin_ph( meas::read_tensor_functions_from_obs<ut::complex>(jObservables(cfg::hedin_ph_imprsum::Worm::name())("susceptibility"), jParams, jParams(cfg::hedin_ph_imprsum::Worm::name()), jHybMatrix, hyb.size()));
            auto const hedin_pp = rearrange_hedin_pp( meas::read_tensor_functions_from_obs<ut::complex>(jObservables(cfg::hedin_pp_imprsum::Worm::name())("susceptibility"), jParams, jParams(cfg::hedin_pp_imprsum::Worm::name()), jHybMatrix, hyb.size()));
            auto const hedin_tph = rearrange_hedin_tph(hedin_ph);
            
            std::vector<io::ctens> vertex = meas::read_tensor_functions_from_obs<ut::complex>(jObservables(cfg::vertex_imprsum::Worm::name())("susceptibility"), jParams, jParams(cfg::vertex_imprsum::Worm::name()), jHybMatrix, hyb.size());
            
            std::cout << "OK" << std::endl;
            
            //Construct interaction matrix U_ijkl
            std::cout << "Constructing interaction matrix ... " << std::endl;
            
            //This matrix has a factor of 1/2 built in, so we must adjust the resulting kernel equations
            params::complete_impurity<Value>(jParams);
            imp::Tensor<Value> const U_tmp(jParams("hloc")("two body"),jHybMatrix.size());
            InteractionTensor<Value> const U(U_tmp,jHybMatrix.size());
            
            std::cout << "OK" << std::endl;
            
            
            if (susc_ph.size() != susc_pp.size())
            throw std::runtime_error("To evaluate asymptotic vertex kernels, all hedin and susc worms must share cutoff frequencies\n");
            
            if (hedin_ph.size() != hedin_pp.size())
            throw std::runtime_error("To evaluate asymptotic vertex kernels, all hedin and susc worms must share cutoff frequencies\n");
            
            if (susc_ph.size() != nMatGB)
            throw std::runtime_error("To evaluate asymptotic vertex kernels, all hedin and susc worms must share cutoff frequencies\n");
            
            
            //Construct Kernel-1 Functions
            std::cout << "Calculating Kernel-1 functions ... " << std::flush;
            
            std::vector<io::ctens> kernel_1_ph(nMatGB_kernel, io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            std::vector<io::ctens> kernel_1_pp(nMatGB_kernel, io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            std::vector<io::ctens> kernel_1_tph(nMatGB_kernel, io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            
            for (int n=0; n<nMatGB_kernel; n++){
                
                int const n_susc_test = omega_b.pos(omega_b_kernel(n));
                bool const do_conj = n_susc_test < 0 ? true : false;
                int const n_susc = do_conj ? omega_b.pos(-omega_b_kernel(n)) : n_susc_test;
                
                for(auto const& ijkl : vertex[0].ijkl()){
                    
                    auto const a=ijkl[1];
                    auto const b=ijkl[2];
                    auto const c=ijkl[3];
                    auto const d=ijkl[4];
                    
                    std::string const entry = std::to_string(2*a)+"_"+std::to_string(2*b+1)+"_"+std::to_string(2*c)+"_"+std::to_string(2*d+1);
                    
                    kernel_1_ph[n].emplace(a,b,c,d,entry,0.);
                    kernel_1_tph[n].emplace(a,b,c,d,entry,0.);
                    kernel_1_pp[n].emplace(a,b,c,d,entry,0.);
                    
                    for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                    for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                    for(std::size_t k = 0; k < jHybMatrix.size(); ++k)
                    for(std::size_t l = 0; l < jHybMatrix.size(); ++l){
                        
                        if (do_conj){
                            kernel_1_ph[n](a,b,c,d) -= 4.*U(a,j,b,i) * std::conj(susc_ph[n_susc].at(i,j,k,l)) * U(l,c,k,d);
                            kernel_1_tph[n](a,b,c,d) -= 4.*U(a,l,i,d) * std::conj(susc_tph[n_susc].at(i,j,k,l)) * U(j,c,b,k);
                            kernel_1_pp[n](a,b,c,d) -= U(a,c,k,i) * std::conj(susc_pp[n_susc].at(i,j,k,l)) * U(l,j,b,d);
                        } else {
                            kernel_1_ph[n](a,b,c,d) -= 4.*U(a,j,b,i) * susc_ph[n_susc].at(i,j,k,l) * U(l,c,k,d);
                            kernel_1_tph[n](a,b,c,d) -= 4.*U(a,l,i,d) * susc_tph[n_susc].at(i,j,k,l) * U(j,c,b,k);
                            kernel_1_pp[n](a,b,c,d) -= U(a,c,k,i) * susc_pp[n_susc].at(i,j,k,l) * U(l,j,b,d);
                        }
                        
                    }
                }
            }
            
            std::cout << "OK" << std::endl;
            
            //Construct Kernel-2 Functions
            std::cout << "Calculating Kernel-2 functions ... " << std::flush;
            
            std::vector<io::ctens> kernel_2_ph(nMatGB_kernel*nMatGF_kernel, io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            std::vector<io::ctens> kernel_2_pp(nMatGB_kernel*nMatGF_kernel, io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            std::vector<io::ctens> kernel_2_tph(nMatGB_kernel*nMatGF_kernel, io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            
            for (int om=0; om<nMatGB_kernel; om++)
            for (int nu=0; nu<nMatGF_kernel; nu++){
                
                int const n = nu + om*nMatGF;
                int const wf = greenOM.pos(omega_f_kernel(nu));
                int const wb = greenOM.pos(omega_f_kernel(nu)-omega_b_kernel(om));
                
                int const n_susc_test = omega_b.pos(omega_b_kernel(om));
                bool const do_conj = n_susc_test < 0 ? true : false;
                int const om_susc = do_conj ? omega_b.pos(-omega_b_kernel(om)) : n_susc_test;
                int const nu_neg = omega_f_kernel.pos(-omega_f_kernel(nu));
                int const n_susc = do_conj ? nu_neg + om_susc*nMatGF : nu + om_susc*nMatGF;
                
                
                for(auto const& ijkl : vertex[0].ijkl()){
                
                    auto const a=ijkl[1];
                    auto const b=ijkl[2];
                    auto const c=ijkl[3];
                    auto const d=ijkl[4];
                
                    std::string const entry = std::to_string(2*a)+"_"+std::to_string(2*b+1)+"_"+std::to_string(2*c)+"_"+std::to_string(2*d+1);
                    
                    kernel_2_ph[n].emplace(a,b,c,d,entry, -kernel_1_ph[om].at(a,b,c,d));
                    kernel_2_tph[n].emplace(a,b,c,d,entry, -kernel_1_tph[om].at(a,b,c,d));
                    kernel_2_pp[n].emplace(a,b,c,d,entry, -kernel_1_pp[om].at(a,b,c,d));
                    
                    for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                    for(std::size_t j = 0; j < jHybMatrix.size(); ++j){
                        
                        if (do_conj){
                            kernel_2_ph[n](a,b,c,d)  -= 2.*std::conj(hedin_ph[n_susc].at(a,b,j,i))*U(i,c,j,d)/(green[wf](a,a)*green[wb](b,b));
                            kernel_2_tph[n](a,b,c,d) -= 2.*std::conj(hedin_tph[n_susc].at(a,i,j,d))*U(i,c,b,j)/(green[wf](a,a)*green[wb](d,d));
                            kernel_2_pp[n](a,b,c,d)  -=    std::conj(hedin_pp[n_susc].at(a,i,c,j))*U(j,i,b,d)/(green[wf](a,a)*green[wb](c,c));
                        } else {
                            kernel_2_ph[n](a,b,c,d)  -= 2.* (hedin_ph[n_susc].at(a,b,j,i))*U(i,c,j,d)/(green[wf](a,a)*green[wb](b,b));
                            kernel_2_tph[n](a,b,c,d) -= 2.* (hedin_tph[n_susc].at(a,i,j,d))*U(i,c,b,j)/(green[wf](a,a)*green[wb](d,d));
                            kernel_2_pp[n](a,b,c,d)  -=     (hedin_pp[n_susc].at(a,i,c,j))*U(j,i,b,d)/(green[wf](a,a)*green[wb](c,c));
                        }
                     
                    }
                }
            }
      
            std::cout << "OK" << std::endl;
            
            //Write results
            std::cout << "Outputting results ... " << std::flush;
            
            
            jsx::value jObservablesOut;
            jsx::value jKernel_1,jKernel_2;
            
            jKernel_1["ph"] = func::write_functions<Value>(jParams, jHybMatrix, kernel_1_ph);
            jKernel_1["tph"] = func::write_functions<Value>(jParams, jHybMatrix, kernel_1_tph);
            jKernel_1["pp"] = func::write_functions<Value>(jParams, jHybMatrix, kernel_1_pp);
            
            jKernel_2["ph"] = func::write_functions<Value>(jParams, jHybMatrix, kernel_2_ph);
            jKernel_2["tph"] = func::write_functions<Value>(jParams, jHybMatrix, kernel_2_tph);
            jKernel_2["pp"] = func::write_functions<Value>(jParams, jHybMatrix, kernel_2_pp);
            
            jObservablesOut["kernel 1"] = std::move(jKernel_1);
            jObservablesOut["kernel 2"] = std::move(jKernel_2);
            
            std::cout << "OK" << std::endl;
            
            return jObservablesOut;
            
        }
        
        
        template<typename Value>
        jsx::value evaluateFullVertexFromKernels(jsx::value jParams, jsx::value const& jObservables) {
            
            auto const name = worm::Kernels::name;
            jsx::value jWorm = jParams(name);
            
            //Full vertex only for +/- range if complex
            int const pos_and_neg_freq = std::is_same<Value,double>::value ? 1 : 2;
            int const nMatGB = pos_and_neg_freq*jWorm("boson cutoff").int64()-pos_and_neg_freq+1;
            int const nMatGF = 2*jWorm("fermion cutoff").int64();
            int const asymptotic_cutoff_l = jParams(name).is("asymptotic cutoff") ? jParams(name)("asymptotic cutoff").int64() : 10;
            int const asymptotic_cutoff_l_4 = std::pow(asymptotic_cutoff_l,4);
            
            func::OmegaMap omega_f(nMatGF,false,true);
            func::OmegaMap omega_b(nMatGB,true,!std::is_same<Value,double>::value);
            
            //Kernels are always for full +/- frequency range
            auto const name_kernel = cfg::hedin_ph_imprsum::Worm::name();
            int const nMatGF_kernel = 2*(jParams(name_kernel)("basis").string() == "matsubara" ? jParams(name_kernel)("fermion cutoff").int64() : ( jParams(name_kernel).is("fermion cutoff") ? jParams(name_kernel)("fermion cutoff").int64() : 50 ));
            int const nMatGB_kernel = 2*jParams(name_kernel)("boson cutoff").int64()-1;
            
            func::OmegaMap omega_f_kernel(nMatGF_kernel,false,true);
            func::OmegaMap omega_b_kernel(nMatGB_kernel,true,true);
            
            jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
            std::vector<io::cmat> hyb; std::vector<io::Matrix<Value>> hybMoments;
            std::tie(hyb, hybMoments) = partition::func::get_hybridisation<Value>(jParams);
            
            params::complete_impurity<Value>(jParams);
            imp::Tensor<Value> const U_tmp(jParams("hloc")("two body"),jHybMatrix.size());
            InteractionTensor<Value> const U(U_tmp,jHybMatrix.size());
            
            //Read in green functions and susc and hedin susceptibilities
            std::cout << "Reading in kernels ... " << std::flush;
            
            std::vector<io::ctens> kernel_1_ph = meas::read_tensor_functions_from_obs<ut::complex>(jObservables(worm::Kernels::name)("kernel 1")("ph"), jParams, jParams(cfg::susc_ph::Worm::name()), jHybMatrix, hyb.size());
            std::vector<io::ctens> kernel_1_tph = meas::read_tensor_functions_from_obs<ut::complex>(jObservables(worm::Kernels::name)("kernel 1")("tph"), jParams, jParams(cfg::susc_ph::Worm::name()), jHybMatrix, hyb.size());
            std::vector<io::ctens> kernel_1_pp = meas::read_tensor_functions_from_obs<ut::complex>(jObservables(worm::Kernels::name)("kernel 1")("pp"), jParams, jParams(cfg::susc_ph::Worm::name()), jHybMatrix, hyb.size());
            
            std::vector<io::ctens> kernel_2_ph = meas::read_tensor_functions_from_obs<ut::complex>(jObservables(worm::Kernels::name)("kernel 2")("ph"), jParams, jParams(cfg::hedin_ph_imprsum::Worm::name()), jHybMatrix, hyb.size());
            std::vector<io::ctens> kernel_2_tph = meas::read_tensor_functions_from_obs<ut::complex>(jObservables(worm::Kernels::name)("kernel 2")("tph"), jParams, jParams(cfg::hedin_ph_imprsum::Worm::name()), jHybMatrix, hyb.size());
            std::vector<io::ctens> kernel_2_pp = meas::read_tensor_functions_from_obs<ut::complex>(jObservables(worm::Kernels::name)("kernel 2")("pp"), jParams, jParams(cfg::hedin_pp_imprsum::Worm::name()), jHybMatrix, hyb.size());
            
            std::cout << "OK" << std::endl;
            
            std::cout << "Reading in measured vertex ... " << std::flush;
            
            std::vector<io::ctens> measured_vertex = meas::read_tensor_functions_from_obs<ut::complex>(jObservables(cfg::vertex_imprsum::Worm::name())("full vertex"), jParams, jParams(cfg::vertex_imprsum::Worm::name()), jHybMatrix, hyb.size());
            
            std::cout << "OK" << std::endl;
            
            //Construct Kernel-2 Functions
            std::cout << "Calculating Asymptotic Vertex from kernels ... " << std::flush;
            
            std::vector<io::ctens> asymptotic_vertex(nMatGB*nMatGF*nMatGF, io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            
            std::cout << "\n";
            
            for (int om=0; om<nMatGB; om++)
            for (int nu1=0; nu1<nMatGF; nu1++)
            for (int nu2=0; nu2<nMatGF; nu2++){
                
                int const nu1_ph = omega_f_kernel.pos(omega_f(nu1));
                int const nu1_tph = nu1_ph;
                int const nu1_pp = nu1_ph;
                
                int const nu2_ph = omega_f_kernel.pos(omega_f(nu2));
                int const nu2_tph = omega_f_kernel.pos(omega_f(nu1) - omega_b(om));
                int const nu2_pp = nu2_ph;
                
                int const om_ph = omega_b_kernel.pos(omega_b(om));
                int const om_tph = omega_b_kernel.pos(omega_f(nu1) - omega_f(nu2));
                int const om_pp = omega_b_kernel.pos(omega_f(nu1) + omega_f(nu2) - omega_b(om));
                
                int const n = nu2 + nu1*nMatGF + om*nMatGF*nMatGF;
                
                int const n1_ph = nu1_ph + om_ph * nMatGF_kernel;
                int const n1_tph = nu1_tph + om_tph * nMatGF_kernel;
                int const n1_pp = nu1_pp + om_pp * nMatGF_kernel;
                
                int const n2_ph = nu2_ph + om_ph * nMatGF_kernel;
                int const n2_tph = nu2_tph + om_tph * nMatGF_kernel;
                int const n2_pp = nu2_pp + om_pp * nMatGF_kernel;
                                
                for(auto const& ijkl : measured_vertex[0].ijkl()){
                
                    auto const a=ijkl[1];
                    auto const b=ijkl[2];
                    auto const c=ijkl[3];
                    auto const d=ijkl[4];
                    
                    asymptotic_vertex[n].emplace(a,b,c,d,measured_vertex[0].entry(a,b,c,d), -2.*U(a,c,b,d));
                    
                    //PH
                    if(1)
                    if (om_ph>=0){
                        asymptotic_vertex[n](a,b,c,d) += kernel_1_ph[om_ph].at(a,b,c,d);
                        if (nu1_tph>=0 and nu2_ph>=0 and 1){
                            //1;
                            asymptotic_vertex[n](a,b,c,d) += kernel_2_ph[n1_ph].at(a,b,c,d) + kernel_2_ph[n2_ph].at(a,b,c,d);// + kernel_1_ph[om_ph](a,b,c,d);
                            
                            
                        } else if (0) {
                            if (nu1_ph>=0)
                                asymptotic_vertex[n](a,b,c,d) += kernel_2_ph[n1_ph].at(a,b,c,d) + 0.5*kernel_1_ph[om_ph].at(a,b,c,d);
                            else if (nu2_ph>=0)
                                asymptotic_vertex[n](a,b,c,d) += kernel_2_ph[n2_ph].at(a,b,c,d) + 0.5*kernel_1_ph[om_ph].at(a,b,c,d);
                            //else
                            //    full_vertex[n](a,b,c,d) -= 2.*kernel_1_ph[om_ph](a,b,c,d);
                        }
                    }
                    
                    //TPH
                    if(1)
                    if (om_tph>=0){
                        asymptotic_vertex[n](a,b,c,d) += kernel_1_tph[om_tph].at(a,b,c,d);
                        if (nu1_tph>=0 and nu2_tph>=0 and 1){
                            
                            asymptotic_vertex[n](a,b,c,d) += kernel_2_tph[n1_tph].at(a,b,c,d) + kernel_2_tph[n2_tph].at(a,b,c,d);
                            
                            
                        } else if (0) {
                            if (nu1_tph>=0)
                                asymptotic_vertex[n](a,b,c,d) += kernel_2_tph[n1_tph].at(a,b,c,d) - 0.5*kernel_1_tph[om_tph].at(a,b,c,d);
                            else if (nu2_tph>=0)
                                asymptotic_vertex[n](a,b,c,d) += kernel_2_tph[n2_tph].at(a,b,c,d) - 0.5*kernel_1_tph[om_tph].at(a,b,c,d);
                            else
                                asymptotic_vertex[n](a,b,c,d) -= 2.*kernel_1_tph[om_tph].at(a,b,c,d);
                        }
                        
                    }
                    
                    //PP
                    if(1)
                    if (om_pp>=0){
                        asymptotic_vertex[n](a,b,c,d) += kernel_1_pp[om_pp].at(a,b,c,d);
                        if (nu1_pp>=0 and nu2_pp>=0 and 1){
                            //1;
                            asymptotic_vertex[n](a,b,c,d) += kernel_2_pp[n1_pp].at(a,b,c,d) + kernel_2_pp[n2_pp].at(a,b,c,d);// + kernel_1_pp[om_pp].at(a,b,c,d);
                            
                            
                        } else if (0){
                            if (nu1_pp>=0)
                                asymptotic_vertex[n](a,b,c,d) += kernel_2_pp[n1_pp].at(a,b,c,d) - kernel_1_pp[om_pp].at(a,b,c,d);
                            else if (nu2_pp>=0)
                                asymptotic_vertex[n](a,b,c,d) += kernel_2_pp[n2_pp].at(a,b,c,d) - kernel_1_pp[om_pp].at(a,b,c,d);
                            else
                                asymptotic_vertex[n](a,b,c,d) -= 3.*kernel_1_pp[om_pp].at(a,b,c,d);
                        }
                    }
                    
                }
            
            }
            
            std::cout << "OK" << std::endl;
            
            std::cout << "Combining asymptotic and measured vertices ... " << std::flush;
            
            func::Frequencies<Value> frequencies_meas(jParams(cfg::vertex_imprsum::Worm::name()));
            auto const& omega_f_meas = frequencies_meas.omega_f();
            auto const& omega_b_meas = frequencies_meas.omega_b();
            
            std::vector<io::ctens> combined_vertex(asymptotic_vertex.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));

            for (int om=0; om<nMatGB; om++)
            for (int nu1=0; nu1<nMatGF; nu1++)
            for (int nu2=0; nu2<nMatGF; nu2++){
                
                int const nu_prod = std::abs(omega_f(nu1) * (omega_f(nu1) - omega_b(om)) * (omega_f(nu2) - omega_b(om)) * omega_f(nu2));
                int const is_static = 1;//!omega_b(om);
                int const is_diag = nu1 == nu2;
                int const use_asymptotic = nu_prod > asymptotic_cutoff_l_4*(is_static + is_diag - is_static*is_diag);
                
                int const om_meas = omega_b_meas.pos(omega_b(om));
                int const nu1_meas = omega_f_meas.pos(omega_f(nu1));
                int const nu2_meas = omega_f_meas.pos(omega_f(nu2));
                
                int const n = nu2 + nu1*nMatGF + om*nMatGF*nMatGF;
                int const n_meas = nu2_meas + nu1_meas*frequencies_meas.nMatGF() + om_meas*frequencies_meas.nMatGF()*frequencies_meas.nMatGF();
                
                for(auto const& ijkl : measured_vertex[0].ijkl()){
                
                    auto const a=ijkl[1];
                    auto const b=ijkl[2];
                    auto const c=ijkl[3];
                    auto const d=ijkl[4];
                    
                    combined_vertex[n].emplace(a,b,c,d, measured_vertex[0].entry(a,b,c,d), 0.);
                    
                    if (om_meas < 0 or nu1_meas < 0 or nu2_meas < 0 or use_asymptotic){
                        combined_vertex[n](a,b,c,d) = asymptotic_vertex[n](a,b,c,d);
                    } else {
                        combined_vertex[n](a,b,c,d) = measured_vertex[n_meas](a,b,c,d);
                    }
                        
                    
                }
                
            }
            
            std::cout << "OK" << std::endl;
            
            
            std::cout << "Outputting results ... " << std::flush;
            
            jsx::value jObservablesOut;
            
            jObservablesOut["full vertex"] = func::write_functions<Value>(jParams, jHybMatrix, combined_vertex);
            jObservablesOut["full vertex (asymptotic)"] = func::write_functions<Value>(jParams, jHybMatrix, asymptotic_vertex);
            
            std::cout << "OK" << std::endl;
            
            return jObservablesOut;
            
        }
        
        
    }
    
}


#endif









