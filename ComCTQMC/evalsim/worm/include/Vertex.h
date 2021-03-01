#ifndef EVALSIM_WORM_INCLUDE_VERTEX
#define EVALSIM_WORM_INCLUDE_VERTEX


#include "../../../include/JsonX.h"
#include "../../../ctqmc/include/config/Worms.h"
#include "functions/Functions.h"
#include "functions/Measurements.h"
#include "functions/Utilities.h"
#include "Common.h"


namespace evalsim {
    
    namespace worm {
        
        namespace func {
        
            namespace vertex {
                
                template <typename Value>
                void compute_disconnected(jsx::value const& jParams, Frequencies<Value> const& frequencies, std::vector<io::cmat> const& green, OmegaMap const& green_OM, std::vector<io::ctens>& disconnected){
                    
                    double const beta = jParams("beta").real64();
                    jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                    
                    auto const nMatGB = frequencies.nMatGB();
                    auto const nMatGF = frequencies.nMatGF();
                    auto const omega_b = frequencies.omega_b();
                    auto const omega_f = frequencies.omega_f();
                    
                    for(std::size_t nu = 0; nu < nMatGB; ++nu)
                        for(std::size_t i_w = 0; i_w < nMatGF; ++i_w)
                            for(std::size_t j_w = 0; j_w < nMatGF; ++j_w){
                                int const n = nu*nMatGF*nMatGF + j_w*nMatGF + i_w;
                                
                                int i_w_g = green_OM.pos(omega_f(i_w));
                                int w = green_OM.pos(omega_f(j_w) - omega_b(nu));
                                
                                if (i_w_g < 0 or w < 0) throw std::runtime_error("vertex:: insufficient frequencies in green's function to compute disconnected part");
                                
                                for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                                    for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                                        for(std::size_t k = 0; k < jHybMatrix.size(); ++k)
                                            for(std::size_t l = 0; l < jHybMatrix.size(); ++l){
                                                
                                               auto const disc = -(
                                                                   ((i == j and k==l and !omega_b(nu)) ? green[i_w_g](i,j)*green[w](k,l) : 0) -
                                                                   ((i == l and j==k and i_w==j_w) ? green[i_w_g](i,l)*green[w](j,k) : 0 )
                                                                   )*beta;

                                                
                                                disconnected[n].emplace(i,j,k,l,
                                                                           "do not use disconnected entries",
                                                                           disc);
                                            }
                            }
                }
                
                
                template<typename Value>
                void compute_and_subtract_disconnected(jsx::value const& jParams, Frequencies<Value> const& frequencies, std::vector<io::cmat> const& green, OmegaMap const& green_OM, std::vector<io::ctens>& full_in_connected_out){
                    
                    double const beta = jParams("beta").real64();
                    jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                    
                    auto const nMatGB = frequencies.nMatGB();
                    auto const nMatGF = frequencies.nMatGF();
                    auto const omega_b = frequencies.omega_b();
                    auto const omega_f = frequencies.omega_f();
                    
                    
                    for(std::size_t nu = 0; nu < nMatGB; ++nu)
                        for(std::size_t i_w = 0; i_w < nMatGF; ++i_w)
                            for(std::size_t j_w = 0; j_w < nMatGF; ++j_w){
                                int const n = nu*nMatGF*nMatGF + j_w*nMatGF + i_w;
                                
                                int i_w_g = green_OM.pos(omega_f(i_w));
                                int w = green_OM.pos(omega_f(j_w)-omega_b(nu));
                                
                                if (i_w_g < 0 or w < 0) throw std::runtime_error("vertex:: insufficient frequencies in green's function to compute disconnected part");
                                
                                for(auto const ijkl : full_in_connected_out[n].ijkl()){
                                
                                auto const i = ijkl[1];
                                auto const j = ijkl[2];
                                auto const k = ijkl[3];
                                auto const l = ijkl[4];
                                    
                                auto const disconnected = -(
                                                            ((i == j and k==l and !omega_b(nu)) ? green[i_w_g](i,j)*green[w](k,l) : 0) -
                                                            ((i == l and j==k and i_w==j_w) ? green[i_w_g](i,l)*green[w](j,k) : 0 )
                                                            )*beta;
                                                
                                //Due to the read function, full is the wrong sign leading to the flipped signs below.
                                full_in_connected_out[n](i,j,k,l) = disconnected - full_in_connected_out[n](i,j,k,l);
                                                
                                }
                            }
                }
                
                
                template<typename Value>
                void compute_connected_from_improved_estimator(jsx::value const& jParams, Frequencies<Value> const& frequencies,
                                                               std::vector<io::cmat> const& green,  std::vector<io::cmat> const& self, OmegaMap const& green_OM,
                                                               std::vector<io::ctens> const& improved_estimator, std::vector<io::ctens>& connected){
                    
                    jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                    
                    std::vector<io::ctens> disconnected(improved_estimator.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
                    compute_disconnected<Value>(jParams,frequencies,green,green_OM,disconnected);
                    
                    auto const nMatGB = frequencies.nMatGB();
                    auto const nMatGF = frequencies.nMatGF();
                    auto const omega_f = frequencies.omega_f();
                    
                    for(std::size_t nu = 0; nu < nMatGB; ++nu)
                        for(std::size_t i_w = 0; i_w < nMatGF; ++i_w)
                            for(std::size_t j_w = 0; j_w < nMatGF; ++j_w){
                                int n = nu*nMatGF*nMatGF+j_w*nMatGF+i_w;
                                
                                int w = green_OM.pos(omega_f(i_w));
                                
                                for(auto const ijkl : improved_estimator[n].ijkl()){
                                    
                                    auto const i = ijkl[1];
                                    auto const j = ijkl[2];
                                    auto const k = ijkl[3];
                                    auto const l = ijkl[4];
                                    
                                    connected[n].emplace(i,j,k,l, improved_estimator[n].entry(i,j,k,l), 0.);
                                    
                                    for(std::size_t m = 0; m < jHybMatrix.size(); ++m){
                                        
                                        connected[n](i,j,k,l) += std::abs(green[w](m,i)) ?
                                                ( green[w](m,i) * self[w](m,i) * disconnected[n](i,j,k,l)
                                                 - green[w](m,i) * improved_estimator[n](i,j,k,l) )/
                                                ( (i==m ? 1.0 : 0.0) + green[w](m,i) * self[w](m,i) )
                                                : 0;
                                                
                                    }
                                }
                            }
                            
                }
                
                template<typename Value>
                void enforce_symmetries(jsx::value const& jParams, Frequencies<Value> const& frequencies, std::vector<io::ctens> const& vertex, std::vector<io::ctens>& symm_vertex){
                    
                    jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                    
                    auto const nMatGB = frequencies.nMatGB();
                    auto const nMatGF = frequencies.nMatGF();
                    auto const omega_b = frequencies.omega_b();
                    auto const omega_f = frequencies.omega_f();
                    
                    bool const have_neg_boson = (std::is_same<Value,double>::value ? false : true);
                    
                    for(std::size_t nu = 0; nu < nMatGB; ++nu)
                        for(std::size_t i_w = 0; i_w < nMatGF; ++i_w)
                            for(std::size_t j_w = 0; j_w < nMatGF; ++j_w){
                                
                                int const n = nu*nMatGF*nMatGF + j_w*nMatGF + i_w;
                                
                                int ik = -1, jl = -1, both = -1;
                                bool complex_ik = false, complex_jl = false, complex_both = have_neg_boson;
                                
                                int j_m_i_freq = omega_f(j_w) - omega_f(i_w);
                                int j_m_i = omega_b.pos(j_m_i_freq);
                                int i_m_j = omega_b.pos(-j_m_i_freq);
                                
                                int i_m_nu = omega_f.pos(omega_f(i_w) - omega_b(nu));
                                int j_m_nu = omega_f.pos(omega_f(j_w) - omega_b(nu));
                                
                                if (have_neg_boson){
                                    
                                    ik = (j_m_nu < 0 or j_m_i < 0) ? -1 : j_m_i*nMatGF*nMatGF + j_w*nMatGF     + j_m_nu;
                                    jl = (i_m_nu < 0 or i_m_j < 0) ? -1 : i_m_j*nMatGF*nMatGF + i_m_nu*nMatGF  + i_w;
                                    both = (i_m_nu < 0 or j_m_nu < 0) ? -1 : nu*nMatGF*nMatGF + i_m_nu*nMatGF  + j_m_nu;
                                    
                                } else {
                                    
                                    int nu_m_j = omega_f.pos(omega_b(nu) - omega_f(j_w));
                                    int nu_m_i = omega_f.pos(omega_b(nu) - omega_f(i_w));
                                    
                                    int const n_nu = omega_b.pos(-omega_b(nu));
                                    
                                    if (!j_m_i_freq){
                                        
                                        ik = (j_m_i < 0 or j_m_nu < 0) ? -1 : j_m_i*nMatGF*nMatGF + j_w*nMatGF  + j_m_nu;
                                        jl = (i_m_j < 0 or i_m_nu < 0) ? -1 : i_m_j*nMatGF*nMatGF + i_m_nu*nMatGF + i_w;
                                        
                                    }
                                    else if (j_m_i_freq < 0){
                                        
                                        complex_ik = true;
                                        
                                        int const n_j_w = omega_f.pos(-omega_f(j_w));
                                        
                                        ik = (i_m_j < 0 or nu_m_j < 0) ? -1 : i_m_j*nMatGF*nMatGF + n_j_w*nMatGF + nu_m_j;
                                        jl = (i_m_j < 0 or i_m_nu < 0) ? -1 : i_m_j*nMatGF*nMatGF + i_m_nu*nMatGF + i_w;
                                        
                                    } else {
                                        
                                        complex_jl = true;
                                        
                                        int const n_i_w = omega_f.pos(-omega_f(i_w));
                                        
                                        ik = (j_m_i < 0 or j_m_nu < 0) ? -1 : j_m_i*nMatGF*nMatGF + j_w*nMatGF  + j_m_nu;
                                        jl = (j_m_i < 0 or nu_m_i < 0) ? -1 : j_m_i*nMatGF*nMatGF + nu_m_i*nMatGF  + n_i_w;
                                        
                                    }
                                    
                                    both = (nu_m_i < 0 or nu_m_j < 0) ? -1 :  n_nu*nMatGF*nMatGF   + nu_m_i*nMatGF + nu_m_j;
                                    
                                }
                                
                                for(auto const ijkl : vertex[n].ijkl()){
                                    
                                    auto const i = ijkl[1];
                                    auto const j = ijkl[2];
                                    auto const k = ijkl[3];
                                    auto const l = ijkl[4];
                                                
                                    symm_vertex[n].emplace(i,j,k,l,
                                                           vertex[n].entry(i,j,k,l),
                                                           vertex[n](i,j,k,l));
                                    
                                    auto const g2_ik = ik < 0 ? 0 :
                                    complex_ik ? std::conj(vertex[ik](k,j,i,l)) : vertex[ik](k,j,i,l);
                                    auto const g2_jl = jl < 0 ? 0 :
                                    complex_jl ? std::conj(vertex[jl](i,l,k,j)) : vertex[jl](i,l,k,j);
                                    auto const g2_both = both < 0 ? 0 :
                                    complex_both ? std::conj(vertex[both](k,l,i,j)) : vertex[both](k,l,i,j);
                                    
                                    if (ik < 0 and jl < 0 and both < 0)
                                        symm_vertex[n](i,j,k,l) = vertex[n](i,j,k,l);
                                    
                                    else if (ik < 0 and jl < 0){
                                        symm_vertex[n](i,j,k,l) = 0.5*(vertex[n](i,j,k,l) + g2_both);
                                        
                                    }
                                    else if (ik < 0 and both < 0){
                                        symm_vertex[n](i,j,k,l) = 0.5*(vertex[n](i,j,k,l) - g2_jl);
                                        
                                    }
                                    else if (both < 0 and jl < 0){
                                        symm_vertex[n](i,j,k,l) = 0.5*(vertex[n](i,j,k,l) - g2_ik);
                                        
                                    }
                                    else if (ik < 0){
                                        symm_vertex[n](i,j,k,l) = 1./3.*(vertex[n](i,j,k,l) - g2_jl + g2_both);
                                        
                                    }
                                    else if (jl < 0){
                                        symm_vertex[n](i,j,k,l) = 1./3.*(vertex[n](i,j,k,l) - g2_ik + g2_both);
                                        
                                    }
                                    else if (both < 0){
                                        symm_vertex[n](i,j,k,l) = 1./3.*(vertex[n](i,j,k,l) - g2_ik - g2_jl);
                                        
                                    }
                                    else{
                                        symm_vertex[n](i,j,k,l) = 0.25*(vertex[n](i,j,k,l) - g2_ik - g2_jl + g2_both);
                                        
                                    }
                                     
                                                    
                                }
                                            
                            }
                }
                
                template<typename Value>
                void compute_full_vertex(jsx::value const& jParams, Frequencies<Value> const& frequencies,
                                         std::vector<io::cmat> const& green, OmegaMap const& green_OM,
                                         std::vector<io::ctens> const& susceptibility, std::vector<io::ctens>& full_vertex){
                    
                    jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                    
                    auto const nMatGB = frequencies.nMatGB();
                    auto const nMatGF = frequencies.nMatGF();
                    auto const omega_b = frequencies.omega_b();
                    auto const omega_f = frequencies.omega_f();
                    
                    for(std::size_t nu = 0; nu < nMatGB; ++nu)
                        for(std::size_t i_w = 0; i_w < nMatGF; ++i_w)
                            for(std::size_t j_w = 0; j_w < nMatGF; ++j_w){
                                
                                int const n = nu*nMatGF*nMatGF + j_w*nMatGF + i_w;
                                
                                int i_w_g = green_OM.pos(omega_f(i_w));
                                int j_w_g = green_OM.pos(omega_f(j_w));
                                int i_m_nu = green_OM.pos(omega_f(i_w) - omega_b(nu));
                                int j_m_nu = green_OM.pos(omega_f(j_w) - omega_b(nu));
                                
                                for(auto const ijkl : susceptibility[n].ijkl()){
                                
                                    auto const i = ijkl[1];
                                    auto const j = ijkl[2];
                                    auto const k = ijkl[3];
                                    auto const l = ijkl[4];
                                    
                                    full_vertex[n].emplace(i,j,k,l,
                                                           susceptibility[n].entry(i,j,k,l),
                                                           susceptibility[n](i,j,k,l)/
                                                            (green[i_w_g](i,i)*green[j_w_g](j,j)*green[i_m_nu](i,i)*green[j_m_nu](j,j))
                                                           );
                                                
                                }
                            }
                }
            }
        }
    }
}


#endif









