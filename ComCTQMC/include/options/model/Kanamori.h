#ifndef INCLUDE_OPTIONS_MODEL_KANAMORI_H
#define INCLUDE_OPTIONS_MODEL_KANAMORI_H

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
        
        struct Kanamori {
            Kanamori() = delete;
            Kanamori(Basis const& basis, jsx::value const& jTwoBody) :
            n_(basis.n()),
            U_(jTwoBody("U").real64()),
            J_(jTwoBody("J").real64()),
            Up_(jTwoBody.is("Uprime") ? jTwoBody("Uprime").real64() : U_ - 2.*J_) {
                if(jTwoBody("approximation").string() != "none")
                    throw std::runtime_error("opt::Kanamori: approximation " + jTwoBody("approximation").string() + " not defined.");
                
            }
            ~Kanamori() = default;
            
            double operator()(int m1, int m2, int m3, int m4) const {
                if(m1 == m2 && m3 == m4 && m1 == m3) return U_;
                if(m1 == m4 && m2 == m3 && m1 != m2) return Up_;
                if(m1 == m3 && m2 == m4 && m1 != m2) return J_;
                if(m1 == m2 && m3 == m4 && m1 != m3) return J_;
                return .0;
            };
            
            int n() const {
                return n_;
            };
            
            std::string approximation() const {
                return "none";
            };
            
        private:
            int const n_;
            double const U_, J_, Up_;
        };
        
    }
    
};

#endif






