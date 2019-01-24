#ifndef SLATERCONDON
#define SLATERCONDON

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

#include "../ClebschGordan/ClebschGordan.h"
#include "../Basis/Basis.h"

#include "../../JsonX.h"
#include "../../IO.h"


namespace Interaction {
    
    struct SlaterCondon {
        SlaterCondon() = delete;
        SlaterCondon(int const l, jsx::value const& jTwoBody) :
        l_(l),
        D_(2*l_ + 1),
        F_(l_ + 1),
        gaunt_((l_+ 1)*D_*D_),
        tensor_(D_*D_*D_*D_),
        approximation_(jTwoBody("approximation").string()) {
            for(int k = 0; k <= l_; ++k)
                F_[k] = jTwoBody("F" + std::to_string(2*k)).real64();
            
            for(int k = 0; k <= l_; ++k)
                for(int m = 0; m < D_; ++m)
                    for(int mp = 0; mp < D_; ++mp)
                        gaunt(k, m, mp) = (2*l_ + 1)/(2*(2*k) + 1.)*cg::clebschGordan(l_, 0, l_, 0, 2*k, 0)*cg::clebschGordan(l_, m - l_, l_, l_ - mp, 2*k, m - mp)*((mp - l_) % 2 ? -1 : 1);
            
            for(int m1 = 0; m1 < D_; ++m1)
                for(int m2 = 0; m2 < D_; ++m2)
                    for(int m3 = 0; m3 < D_; ++m3)
                        for(int m4 = 0; m4 < D_; ++m4)
                            tensor_[D_*D_*D_*m1 +
                                    D_*D_*m2 +
                                    D_*m3 +
                                    m4] = tensor(m1, m2, m3, m4);
            
        };
        
        double operator()(int m1, int m2, int m3, int m4) const {
            return tensor_[D_*D_*D_*m1 +
                           D_*D_*m2 +
                           D_*m3 +
                           m4];
        };
        
        std::string approximation() const {
            return approximation_;
        };
        
        ~SlaterCondon() = default;
    private:
        int const l_;
        int const D_;
        std::vector<double> F_;
        std::vector<double> gaunt_; // there may be a global phase difference between these "gaunts" and the ones defined in litterature, did not check. haha, who cares.
        std::vector<double> tensor_;
        
        std::string const approximation_;
        
        double& gaunt(int k, int m, int mp) {
            return gaunt_[k*D_*D_ + m*D_ + mp];
        };
        
        double tensor(int m1, int m2, int m3, int m4) {
            double temp = .0;
            
            if(m1 + m2 == m3 + m4)
                for(int k = 0; k <= l_; ++k)
                    temp += F_[k]*gaunt(k, m1, m4)*gaunt(k, m3, m2);
            
            return temp;
        };
    };
};

#endif






