#ifndef INCLUDE_OPTIONS_MODEL_BASIS_H
#define INCLUDE_OPTIONS_MODEL_BASIS_H

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

#include "../../JsonX.h"


namespace opt {
    
    namespace model {
        
        struct Basis {
            Basis() = delete;
            Basis(jsx::value const& jBasis) : n_(jBasis("orbitals").int64()) {
                io::rvec temp(2*n_, 1.);
                
                qns_["N"] = temp;
                
                for(int i =  0; i   < n_; ++i) temp[i] =  .5;
                for(int i = n_; i < 2*n_; ++i) temp[i] = -.5;
                qns_["Sz"] = temp;
            };
            ~Basis() = default;
            
            std::complex<double> operator()(int ms, int m, int s) const {
                return ( ms%n_ == m )&&( ms/n_ == s ) ? 1. : 0;
            }
            
            int n() const {
                return n_;
            };
            
            std::map<std::string, io::rvec> const& qns() const {
                return qns_;
            };
        private:
            int const n_;
            std::map<std::string, io::rvec> qns_;
        };
        
    }
    
};

#endif






