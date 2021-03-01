#ifndef INCLUDE_OPTIONS_SPHERICALHARMONICS_BASIS_H
#define INCLUDE_OPTIONS_SPHERICALHARMONICS_BASIS_H

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
        
        struct Basis {
            Basis() = delete;
            Basis(jsx::value const& jBasis) {
                std::map<std::string, int> lOption{{"s", 0}, {"p", 1}, {"d", 2}, {"f", 3}};
                
                if(!lOption.count(jBasis("orbitals").string()))
                    throw std::runtime_error("opt: orbitals " + jBasis("orbitals").string() + " not defined");
                
                l_ = lOption.at(jBasis("orbitals").string());
            }
            
            int l() const {
                return l_;
            }
            
            int n() const {
                return 2*l_ + 1;
            }
            
        private:
            int l_;
        };
        
        
        struct Imag : Basis {
            Imag() = delete;
            Imag(jsx::value const& jBasis) : Basis(jBasis) {
                io::rvec temp(4*l() + 2, 1.);
                
                qns_["N"] = temp;
                
                for(int i =         0; i < 2*l() + 1; ++i) temp[i] =  .5;
                for(int i = 2*l() + 1; i < 4*l() + 2; ++i) temp[i] = -.5;
                qns_["Sz"] = temp;
                
                for(int i = 0; i < 2*l() + 1; ++i) temp[i] = temp[i + 2*l() + 1] = i - l();
                qns_["Lz"] = temp;
            }
            ~Imag() = default;
            
            std::complex<double> operator()(int ms, int m, int s) const {
                return ( ms%(2*l() + 1) == m )&&( ms/(2*l() + 1) == s ) ? 1. : 0;
            }
            
            std::map<std::string, io::rvec> const& qns() const {
                return qns_;
            }
            
        private:
            std::map<std::string, io::rvec> qns_;
        };
        
        
        struct Real : Basis {
            Real() = delete;
            Real(jsx::value const& jBasis) : Basis(jBasis) {
                io::rvec temp(4*l() + 2, 1.);
                
                qns_["N"] = temp;
                
                for(int i =         0; i < 2*l() + 1; ++i) temp[i] =  .5;
                for(int i = 2*l() + 1; i < 4*l() + 2; ++i) temp[i] = -.5;
                qns_["Sz"] = temp;
            }
            
            std::complex<double> operator()(int os, int m, int s) const {
                if(os/(2*l() + 1) != s) return .0;
                
                int o = os%(2*l() + 1);
                
                o -= l();
                m -= l();
                
                if(o == 0 && m == 0) return 1.;
                
                std::complex<double> const I(.0, 1.);
                
                if(o < 0) {
                    if(m ==  o) return                I/std::sqrt(2.);
                    if(m == -o) return -(o%2 ? -I  : I)/std::sqrt(2.);
                } else {
                    if(m == -o) return                1./std::sqrt(2.);
                    if(m ==  o) return  (o%2 ? -1. : 1.)/std::sqrt(2.);
                }
                
                return .0;
            }
            
            std::map<std::string, io::rvec> const& qns() const {
                return qns_;
            }
            
        private:
            std::map<std::string, io::rvec> qns_;
        };
        
        
        struct Coupled : Basis {
            Coupled() = delete;
            Coupled(jsx::value const& jBasis) : Basis(jBasis) {
                io::rvec temp(4*l() + 2, 1.);
                
                qns_["N"] = temp;
                
                for(int i =     0; i <      2*l(); ++i) temp[i] = i -   l() + .5;
                for(int i = 2*l(); i <  4*l() + 2; ++i) temp[i] = i - 3*l() - .5;
                qns_["Jz"] = temp;
            }
            ~Coupled() = default;
            
            std::complex<double> operator()(int jm_j, int m, int s) const {
                double j   = jm_j < 2*l() ?        l() - .5 :          l() + .5;
                double m_j = jm_j < 2*l() ? jm_j - l() + .5 : jm_j - 3*l() - .5;
                
                return cg::clebschGordan(l(), m - l(), 0.5, s - 0.5, j, m_j); // m s oder s m ????
            }
            
            std::map<std::string, io::rvec> const& qns() const {
                return qns_;
            }
            
        private:
            std::map<std::string, io::rvec> qns_;
        };
        
    }
    
};

#endif






