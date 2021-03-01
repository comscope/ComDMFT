#ifndef INCLUDE_OPTIONS_MODEL_OBSERVABLES_H
#define INCLUDE_OPTIONS_MODEL_OBSERVABLES_H

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
    
    namespace model {
        
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
        
    }
    
};

#endif






