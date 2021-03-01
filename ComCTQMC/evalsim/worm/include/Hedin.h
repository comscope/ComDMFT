#ifndef EVALSIM_WORM_INCLUDE_HEDIN
#define EVALSIM_WORM_INCLUDE_HEDIN


#include "../../../include/JsonX.h"
#include "../../../ctqmc/include/config/Worms.h"
#include "functions/Functions.h"
#include "functions/Measurements.h"
#include "functions/Utilities.h"
#include "Common.h"


namespace evalsim {
    
    namespace worm {
        
        namespace func {
        
        namespace hedin {
            
            template<typename Value>
            void impl_compute_connected_from_improved_estimator(jsx::value const& jParams, Frequencies<Value> const& frequencies,
                                                                std::vector<io::cmat> const& green,  std::vector<io::cmat> const& self, OmegaMap const& green_OM,
                                                                std::vector<io::ctens> const& improved_estimator, std::vector<io::ctens>const& disconnected, std::vector<io::ctens>& connected){
                
                jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                
                auto const nMatGB = frequencies.nMatGB();
                auto const nMatGF = frequencies.nMatGF();
                auto const omega_f = frequencies.omega_f();
                
                for(std::size_t nu = 0; nu < nMatGB; ++nu)
                    for(std::size_t i_w = 0; i_w < nMatGF; ++i_w){
                        int n = nu*nMatGF + i_w;
                        int w = green_OM.pos(omega_f(i_w));
                        
                        for(auto const ijkl : improved_estimator[n].ijkl()){
                        
                            auto const i = ijkl[1];
                            auto const j = ijkl[2];
                            auto const k = ijkl[3];
                            auto const l = ijkl[4];
                             
                            connected[n].emplace(i,j,k,l,
                                                 improved_estimator[n].entry(i,j,k,l),
                                                 0.0);
                            
                            for(std::size_t m = 0; m < jHybMatrix.size(); ++m){
                                
                                connected[n](i,j,k,l) += std::abs(green[w](m,i)) ?
                                (green[w](m,i) * self[w](m,i) * disconnected[n](i,j,k,l)
                                 - green[w](m,i) * improved_estimator[n](i,j,k,l) )/
                                ( (i==m ? 1.0 : 0.0) + green[w](m,i) * self[w](m,i) )
                                : 0;

                            }
                        }
                    }
            }
            
            
            
            namespace ph {
            
                template <typename Value>
                void compute_disconnected(jsx::value const& jParams, Frequencies<Value> const& frequencies,
                                          std::vector<io::cmat> const& green, OmegaMap const& green_OM,
                                          jsx::value const& jOccupation, std::vector<io::ctens>& disconnected){
                    
                    double const beta = jParams("beta").real64();
                    jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                    
                    auto const nMatGB = frequencies.nMatGB();
                    auto const nMatGF = frequencies.nMatGF();
                    auto const omega_b = frequencies.omega_b();
                    auto const omega_f = frequencies.omega_f();
                    
                    for(std::size_t i_w = 0; i_w < nMatGF; ++i_w)
                        for(std::size_t nu = 0; nu < nMatGB; ++nu){
                            int n = i_w+nu*nMatGF;
                            
                            int i_w_g = green_OM.pos(omega_f(i_w));
                            int w = green_OM.pos(omega_f(i_w)-omega_b(nu));
                            
                            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                                for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                                    for(std::size_t k = 0; k < jHybMatrix.size(); ++k)
                                        for(std::size_t l = 0; l < jHybMatrix.size(); ++l){
                                            
                                            std::string const entry = jHybMatrix(l)(l).string();
                                            auto const& occ = jsx::at<io::Vector<Value>>(jOccupation(entry));
                                            
                                            auto const disc = green[i_w_g](i,i)*( beta*((i==j and k==l and !omega_b(nu)) ? (occ[0]) : 0.0)
                                                                     - ((i==k and l==j) ? green[w](l,l) : 0 ));
                                            
                                            if (disc.real() or disc.imag())
                                                disconnected[n].emplace(i,j,k,l,
                                                                      "do not use disconnected entries",
                                                                      disc);
                                            
                                        }
                        }
                    
                }
                
                template <typename Value>
                void compute_and_subtract_disconnected(jsx::value const& jParams, Frequencies<Value> const& frequencies,
                                                       std::vector<io::cmat> const& green, OmegaMap const& green_OM,
                                                       jsx::value const& jOccupation, std::vector<io::ctens>& full_in_connected_out){
                    
                    jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                    
                    double const beta = jParams("beta").real64();
                    auto const nMatGB = frequencies.nMatGB();
                    auto const nMatGF = frequencies.nMatGF();
                    auto const omega_b = frequencies.omega_b();
                    auto const omega_f = frequencies.omega_f();
                    
                    for(std::size_t i_w = 0; i_w < nMatGF; ++i_w)
                        for(std::size_t nu = 0; nu < nMatGB; ++nu){
                            int n = i_w+nu*nMatGF;
                            int i_w_g = green_OM.pos(omega_f(i_w));
                            int w = green_OM.pos(omega_f(i_w)-omega_b(nu));
                            
                            for(auto const ijkl : full_in_connected_out[n].ijkl()){
                            
                                auto const i = ijkl[1];
                                auto const j = ijkl[2];
                                auto const k = ijkl[3];
                                auto const l = ijkl[4];
                                            
                                full_in_connected_out[n](i,j,k,l)*=-1.0;
                                            
                                std::string const entry = jHybMatrix(l)(l).string();
                                auto const& occ = jsx::at<io::Vector<Value>>(jOccupation(entry));
                                            
                                full_in_connected_out[n](i,j,k,l) += green[i_w_g](i,i)*( beta*((i==j and k==l and !omega_b(nu)) ? (occ[0]) : 0.0)
                                                                                              - ((i==k and l==j) ? green[w](l,l) : 0 ));
                                //susc = disc - full
                            }
                        }
                }
                
                template <typename Value>
                void compute_connected_from_improved_estimator(jsx::value const& jParams, Frequencies<Value> const& frequencies,
                                                               std::vector<io::cmat> const& green,  std::vector<io::cmat> const& self, OmegaMap const& green_OM,
                                                               jsx::value const& jOccupation, std::vector<io::ctens> const& improved_estimator, std::vector<io::ctens>& connected){
                    
                    jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                    
                    std::vector<io::ctens> disconnected(improved_estimator.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
                    compute_disconnected<Value>(jParams, frequencies, green, green_OM, jOccupation, disconnected);
                    
                    impl_compute_connected_from_improved_estimator(jParams, frequencies, green, self, green_OM, improved_estimator, disconnected, connected);
                    
                }
                
                template <typename Value>
                void enforce_symmetries(jsx::value const& jParams, Frequencies<Value> const& frequencies, std::vector<io::ctens> const& no_symm, std::vector<io::ctens>& symm){
                    
                    jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                    
                    auto const nMatGB = frequencies.nMatGB();
                    auto const nMatGF = frequencies.nMatGF();
                    auto const omega_b = frequencies.omega_b();
                    auto const omega_f = frequencies.omega_f();
                    
                    for(std::size_t i_w = 0; i_w < nMatGF; ++i_w)
                        for(std::size_t nu = 0; nu < nMatGB; ++nu){
                            int n = i_w + nu*nMatGF;
                            
                            //swapping operator moves nu -> omega-nu
                            int m = omega_f.pos(omega_b(nu)-omega_f(i_w));
                            if (m >= 0) m += nu*nMatGF;
                            
                            for(auto const ijkl : no_symm[n].ijkl()){
                            
                                auto const i = ijkl[1];
                                auto const j = ijkl[2];
                                auto const k = ijkl[3];
                                auto const l = ijkl[4];
                                            
                                symm[n].emplace(i,j,k,l,
                                                no_symm[n].entry(i,j,k,l),
                                                0.);
                                //Original + hermetian symm
                                if (m >=0){
                                    symm[n](i,j,k,l) = 0.5*(no_symm[n](i,j,k,l) + std::conj(no_symm[m](j,i,l,k)));
                                } else {
                                    symm[n](i,j,k,l) = (no_symm[n](i,j,k,l));
                                }
                                            
                            }
                        }
                }
            }
            
            namespace pp {
                
                template <typename Value>
                void compute_disconnected(jsx::value const& jParams, Frequencies<Value> const& frequencies,
                                          std::vector<io::cmat> const& green, OmegaMap const& green_OM,
                                          jsx::value const& jOccupation, std::vector<io::ctens>& disconnected){
                    
                    jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                    
                    auto const nMatGB = frequencies.nMatGB();
                    auto const nMatGF = frequencies.nMatGF();
                    auto const omega_b = frequencies.omega_b();
                    auto const omega_f = frequencies.omega_f();
                    
                    for(std::size_t i_w = 0; i_w < nMatGF; ++i_w)
                        for(std::size_t nu = 0; nu < nMatGB; ++nu){
                            int n = i_w+nu*nMatGF;
                            int i_w_g = green_OM.pos(omega_f(i_w));
                            int w = green_OM.pos(omega_b(nu)-omega_f(i_w));
                            
                            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                                for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                                    for(std::size_t k = 0; k < jHybMatrix.size(); ++k)
                                        for(std::size_t l = 0; l < jHybMatrix.size(); ++l){
                                            
                                            auto const disc = ((i==k and j==l ? 1. : 0.) - (i==l and j==k ? 1. : 0.))*green[i_w_g](i,i)*green[w](j,j);
                                            
                                            if (disc.real() or disc.imag())
                                                 disconnected[n].emplace(i,j,k,l,
                                                                       "do not use disconnected entries",
                                                                       disc);
                                        }
                        }
                    
                }
                
                template <typename Value>
                void compute_and_subtract_disconnected(jsx::value const& jParams, Frequencies<Value> const& frequencies,
                                                       std::vector<io::cmat> const& green, OmegaMap const& green_OM, jsx::value const& jOccupation,
                                                       std::vector<io::ctens>& full_in_connected_out){
                    
                    jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                    
                    auto const nMatGB = frequencies.nMatGB();
                    auto const nMatGF = frequencies.nMatGF();
                    auto const omega_b = frequencies.omega_b();
                    auto const omega_f = frequencies.omega_f();
                    
                    for(std::size_t i_w = 0; i_w < nMatGF; ++i_w)
                        for(std::size_t nu = 0; nu < nMatGB; ++nu){
                            int n = i_w+nu*nMatGF;
                            int i_w_g = green_OM.pos(omega_f(i_w));
                            int w = green_OM.pos(omega_b(nu)-omega_f(i_w));
                            
                            for(auto const ijkl : full_in_connected_out[n].ijkl()){
                            
                                auto const i = ijkl[1];
                                auto const j = ijkl[2];
                                auto const k = ijkl[3];
                                auto const l = ijkl[4];
                                            
                                full_in_connected_out[n](i,j,k,l)*=-1.0;
                                            
                                full_in_connected_out[n](i,j,k,l) += ((i==k and j==l ? 1. : 0.) - (i==l and j==k ? 1. : 0.))*green[i_w_g](i,i)*green[w](j,j);
                                            
                            }
                        }
                    
                }
                

                template <typename Value>
                void compute_connected_from_improved_estimator(jsx::value const& jParams, Frequencies<Value> const& frequencies,
                                                               std::vector<io::cmat> const& green,  std::vector<io::cmat> const& self,
                                                               OmegaMap const& green_OM, jsx::value const& jOccupation,
                                                               std::vector<io::ctens> const& improved_estimator, std::vector<io::ctens>& connected){
                    
                    jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                    
                    std::vector<io::ctens> disconnected(improved_estimator.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
                    compute_disconnected<Value>(jParams, frequencies, green, green_OM, jOccupation, disconnected);
                    
                    impl_compute_connected_from_improved_estimator(jParams, frequencies, green, self, green_OM, improved_estimator, disconnected, connected);
                    
                }
                
                template <typename Value>
                void enforce_symmetries(jsx::value const& jParams, Frequencies<Value> const& frequencies, std::vector<io::ctens> const& no_symm, std::vector<io::ctens>& symm){
                    
                    jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                    
                    auto const nMatGB = frequencies.nMatGB();
                    auto const nMatGF = frequencies.nMatGF();
                    auto const omega_b = frequencies.omega_b();
                    auto const omega_f = frequencies.omega_f();
                    
                    for(std::size_t i_w = 0; i_w < nMatGF; ++i_w)
                        for(std::size_t nu = 0; nu < nMatGB; ++nu){
                            int n = i_w + nu*nMatGF;
                            
                            //swapping operator moves nu -> omega-nu
                            int m = omega_f.pos(omega_b(nu)-omega_f(i_w));
                            if (m >= 0) m += nu*nMatGF;
                            
                            for(auto const ijkl : no_symm[n].ijkl()){
                            
                                auto const i = ijkl[1];
                                auto const j = ijkl[2];
                                auto const k = ijkl[3];
                                auto const l = ijkl[4];
                                
                                symm[n].emplace(i,j,k,l,
                                                no_symm[n].entry(i,j,k,l),
                                                0.);
                                            
                                //Original + bilinear swap + operator swap + both swaps ()
                                if (m >=0 ){
                                    symm[n](i,j,k,l) = 0.25*(no_symm[n](i,j,k,l) - no_symm[n](i,j,l,k) - no_symm[m](j,i,k,l) + no_symm[m](j,i,l,k));
                                //Original + bilinear swap
                                } else {
                                    symm[n](i,j,k,l) = 0.5*(no_symm[n](i,j,k,l) - no_symm[n](i,j,l,k));
                                }
                                                
                            }
                        }
                }
                
            }
            
            
        }
        
        }
        
    }
    
}


#endif









