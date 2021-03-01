#ifndef INCLUDE_OPTIONS_SPHERICALHARMONICS_SLATERCONDON_H
#define INCLUDE_OPTIONS_SPHERICALHARMONICS_SLATERCONDON_H

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

#include "ClebschGordan.h"

#include "../../JsonX.h"

namespace opt {
    
    namespace sphericalharmonics {
        
        struct SlaterCondon {
            SlaterCondon() = delete;
            SlaterCondon(Basis const& basis, jsx::value const& jTwoBody) :
            l_(basis.l()), n_(basis.n()),
            F_(l_ + 1),
            gaunt_((l_+ 1)*n_*n_),
            approximation_(jTwoBody("approximation").string()) {
                if(approximation_ != "ising" && approximation_ != "none")
                    throw std::runtime_error("opt::SlaterCondon: approximation " + approximation_ + " not defined.");
                
                for(int k = 0; k <= l_; ++k)
                    F_[k] = jTwoBody("F" + std::to_string(2*k)).real64();
                
                for(int k = 0; k <= l_; ++k)
                    for(int m = 0; m < n_; ++m)
                        for(int mp = 0; mp < n_; ++mp)
                            gaunt(k, m, mp) = (2*l_ + 1)/(2*(2*k) + 1.)*cg::clebschGordan(l_, 0, l_, 0, 2*k, 0)*cg::clebschGordan(l_, m - l_, l_, l_ - mp, 2*k, m - mp)*((mp - l_) % 2 ? -1 : 1);
            };
            ~SlaterCondon() = default;
            
            double operator()(int m1, int m2, int m3, int m4) const {
                double temp = .0;
                
                if(m1 + m2 == m3 + m4)
                    for(int k = 0; k <= l_; ++k)
                        temp += F_[k]*gaunt(k, m1, m4)*gaunt(k, m3, m2);
                
                return temp;
            };
            
            int n() const {
                return n_;
            };
            
            std::string const& approximation() const {
                return approximation_;
            };
            
        private:
            int const l_, n_;
            std::vector<double> F_;
            std::vector<double> gaunt_; // there may be a global phase difference between these "gaunts" and the ones defined in litterature, did not check. haha, who cares.
            std::string const approximation_;
            
            double& gaunt(int k, int m, int mp) {
                return gaunt_[k*n_*n_ + m*n_ + mp];
            };
            
            double gaunt(int k, int m, int mp) const {
                return gaunt_[k*n_*n_ + m*n_ + mp];
            };
        };
        
    }
    
};

#endif






