#ifndef EVALSIM_INCLUDE_FUNCTIONS_UTILITIES_H
#define EVALSIM_INCLUDE_FUNCTIONS_UTILITIES_H

#include <tuple>
#include <vector>
#include <string>

#include "../../../../include/measurements/Measurements.h"
#include "../../../../include/io/Matrix.h"
#include "../../../../ctqmc/include/Utilities.h"


namespace evalsim {
    
    namespace worm {
        
        namespace func {
            
            struct iOmega {
                iOmega() = delete;
                iOmega(double const beta, int const is_fermionic=1) : beta_(beta), is_fermionic_(is_fermionic) {};
                std::complex<double> operator()(int n) const {
                    return {.0, (2*n + is_fermionic_)*M_PI/beta_};
                };
            private:
                double const beta_;
                int const is_fermionic_;
            };
            
        
            struct OmegaMap{
                
                OmegaMap() = delete;
                OmegaMap(std::size_t const n, bool bosonic, bool symmetric) : n_(n), val_(n,0){
                    
                    if (symmetric){
                        if (!bosonic and !n%2) throw std::runtime_error("symmetric OmegaMap: number of fermionic frequencies must be even\n");
                        if (bosonic and n%2 == n) throw std::runtime_error("symmetric OmegaMap: number of bosonic frequencies must be odd\n");
                    }
                    
                    int const start = !symmetric ? 0 : -n_/2;
                    int const shift = bosonic ? 0 : 1;
                    
                    for (std::size_t i=0;i<n_; i++){
                        val_[i]=2*(start+i)+shift;
                        pos_[val_[i]]=i;
                    }
                    
                }
                
                int operator()(std::size_t const i) const { return val_[i];}
                inline std::size_t pos(int const i) const { auto it = pos_.find(i); return pos_.end() == it ? -1 : it->second;}
                
            private:
                
                int const n_;
                
                std::map<int,std::size_t> pos_;
                std::vector<int> val_;
                
            };
        
            std::vector<io::cmat> green_function_on_full_axis(std::vector<io::cmat>const& green){
                std::vector<io::cmat> r(2*green.size(),io::cmat(green[0].I(),green[0].J()));
                
                auto it_forward = r.begin() + green.size();
                auto it_backward = r.begin() + green.size()-1;
                for(auto const& x : green){
                    *it_forward++ = x;
                    *it_backward-- = x.conj();
                }
                    
                return r;
            }
            
        }
        
    }
    
}

#endif //EVALSIM










