#ifndef OPT_INTERACTION_H
#define OPT_INTERACTION_H

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

#include "../JsonX.h"


namespace opt {
    
    struct SlaterCondon {
        SlaterCondon() = delete;
        SlaterCondon(int const l, jsx::value const& jTwoBody) :
        l_(l),
        n_(2*l_ + 1),
        F_(l_ + 1),
        gaunt_((l_+ 1)*n_*n_),
        tensor_(n_*n_*n_*n_),
        approximation_(jTwoBody("approximation").string()){
            if(approximation_ != "ising" && approximation_ != "none")
                throw std::runtime_error("opt::SlaterCondon: approximation " + approximation_ + " not defined.");
            
            for(int k = 0; k <= l_; ++k)
                F_[k] = jTwoBody("F" + std::to_string(2*k)).real64();
            
            for(int k = 0; k <= l_; ++k)
                for(int m = 0; m < n_; ++m)
                    for(int mp = 0; mp < n_; ++mp)
                        gaunt(k, m, mp) = (2*l_ + 1)/(2*(2*k) + 1.)*cg::clebschGordan(l_, 0, l_, 0, 2*k, 0)*cg::clebschGordan(l_, m - l_, l_, l_ - mp, 2*k, m - mp)*((mp - l_) % 2 ? -1 : 1);
            
            for(int m1 = 0; m1 < n_; ++m1)
                for(int m2 = 0; m2 < n_; ++m2)
                    for(int m3 = 0; m3 < n_; ++m3)
                        for(int m4 = 0; m4 < n_; ++m4)
                            tensor_[n_*n_*n_*m1 +
                                    n_*n_*m2 +
                                    n_*m3 +
                                    m4] = tensor(m1, m2, m3, m4);
        };
        ~SlaterCondon() = default;
        
        double operator()(int m1, int m2, int m3, int m4) const {
            return tensor_[n_*n_*n_*m1 +
                           n_*n_*m2 +
                           n_*m3 +
                           m4];
        };
        
        int n() const { return n_;};
        std::string const& approximation() const { return approximation_;};
    private:
        int const l_, n_;
        std::vector<double> F_;
        std::vector<double> gaunt_; // there may be a global phase difference between these "gaunts" and the ones defined in litterature, did not check. haha, who cares.
        std::vector<double> tensor_;
        std::string const approximation_;
        
        double& gaunt(int k, int m, int mp) {
            return gaunt_[k*n_*n_ + m*n_ + mp];
        };
        
        double tensor(int m1, int m2, int m3, int m4) {
            double temp = .0;
            
            if(m1 + m2 == m3 + m4)
                for(int k = 0; k <= l_; ++k)
                    temp += F_[k]*gaunt(k, m1, m4)*gaunt(k, m3, m2);
            
            return temp;
        };
    };
    
    
    
    struct Kanamori {
        Kanamori() = delete;
        Kanamori(int n, jsx::value const& jTwoBody) :
        n_(n),
        U_(jTwoBody("U").real64()),
        J_(jTwoBody("J").real64()),
        Up_(jTwoBody.is("Uprime") ? jTwoBody("Uprime").real64() : U_ - 2.*J_),
        tensor_(n_*n_*n_*n_) {
            if(jTwoBody("approximation").string() != "none")
                throw std::runtime_error("opt::Kanamori: approximation " + jTwoBody("approximation").string() + " not defined.");
            
            for(int m1 = 0; m1 < n_; ++m1)
                for(int m2 = 0; m2 < n_; ++m2)
                    for(int m3 = 0; m3 < n_; ++m3)
                        for(int m4 = 0; m4 < n_; ++m4)
                            tensor_[n_*n_*n_*m1 +
                                    n_*n_*m2 +
                                    n_*m3 +
                                    m4] = tensor(m1, m2, m3, m4);
        }
        ~Kanamori() = default;
        
        double operator()(int m1, int m2, int m3, int m4) const {
            return tensor_[n_*n_*n_*m1 +
                           n_*n_*m2 +
                           n_*m3 +
                           m4];
        };
        
        int n() const { return n_;};
        std::string approximation() const { return "none";};
    private:
        int const n_;
        double const U_, J_, Up_;
        std::vector<double> tensor_;
        
        double tensor(int m1, int m2, int m3, int m4) const {
            if(m1 == m2 && m3 == m4 && m1 == m3) return U_;
            if(m1 == m4 && m2 == m3 && m1 != m2) return Up_;
            if(m1 == m3 && m2 == m4 && m1 != m2) return J_;
            if(m1 == m2 && m3 == m4 && m1 != m3) return J_;
            return .0;
        };
    };
};

#endif






