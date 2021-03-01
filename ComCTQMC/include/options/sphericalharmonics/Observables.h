#ifndef INCLUDE_OPTIONS_SPHERICALHARMONICS_OBSERVABLES_H
#define INCLUDE_OPTIONS_SPHERICALHARMONICS_OBSERVABLES_H

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <bitset>
#include <cassert>
#include <iomanip>
#include <set>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <string>

#include "Basis.h"

namespace opt {
    
    namespace sphericalharmonics {
        
        struct S2 {
            S2() = delete;
            S2(Basis const& basis) : n_(basis.n()) {};
            ~S2() = default;
            
            int N() const { return 2*n_;};
            
            double t(int fDagg, int f) const {
                return .0;
            };
            double V(int f1Dagg, int f1, int f2Dagg, int f2) const {
                double temp = .0;
                if(f1Dagg%n_ == f1%n_ && f2Dagg%n_ == f2%n_) {
                    if(f1Dagg/n_ == f2/n_ && f1/n_ == f2Dagg/n_) temp += 2.;
                    if(f1Dagg/n_ == f1/n_ && f2Dagg/n_ == f2/n_) temp -= 1.;
                }
                return temp/4.;
            }
            
        private:
            int const n_;
        };
        
        struct J2 {
            J2() = delete;
            J2(Basis const& basis) : l_(basis.l()), Jz_(4*l_ + 2, 4*l_ + 2), Jp_(4*l_ + 2, 4*l_ + 2) {
                for(int i =    0; i <     2*l_; ++i) Jz_(i, i) = i - l_ + .5;
                for(int i = 2*l_; i < 4*l_ + 2; ++i) Jz_(i, i) = i - 3*l_ - .5;
                
                for(int i =    0; i < 2*l_ - 1; ++i) Jp_(i + 1, i) = std::sqrt((l_ - .5)*(l_ +  .5) - (i - l_ + .5)*(i - l_ + .5 + 1.));
                for(int i = 2*l_; i < 4*l_ + 1; ++i) Jp_(i + 1, i) = std::sqrt((l_ + .5)*(l_ + 1.5) - (i - 3*l_ - .5)*(i - 3*l_ - .5 + 1.));
            };
            ~J2() = default;
            
            int N() const { return 4*l_ + 2;};
            
            double t(int fDagg, int f) const {
                return .0;
            };
            
            double V(int f1Dagg, int f1, int f2Dagg, int f2) const {
                return Jz_(f1Dagg, f1)*Jz_(f2Dagg, f2) + 0.5*Jp_(f1Dagg, f1)*Jp_(f2, f2Dagg) + 0.5*Jp_(f1, f1Dagg)*Jp_(f2Dagg, f2);
            }
            
        private:
            int const l_;
            io::rmat Jz_;
            io::rmat Jp_;
        };
        
    }
    
};

#endif






