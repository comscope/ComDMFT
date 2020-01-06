#ifndef OPT_BASIS_H
#define OPT_BASIS_H

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
    
    struct SphericalImag {
        SphericalImag() = delete;
        SphericalImag(int l) : l_(l) {
            io::rvec temp(4*l_ + 2, 1.);
            
            qns_["N"] = temp;
            
            for(int i =       0; i < 2*l + 1; ++i) temp[i] =  .5;
            for(int i = 2*l + 1; i < 4*l + 2; ++i) temp[i] = -.5;
            qns_["Sz"] = temp;
            
            for(int i = 0; i < 2*l + 1; ++i) temp[i] = temp[i + 2*l + 1] = i - l_;
            qns_["Lz"] = temp;
        }
        ~SphericalImag() = default;
        
        std::complex<double> operator()(int ms, int m, int s) const {
            return ( ms%(2*l_ + 1) == m )&&( ms/(2*l_ + 1) == s ) ? 1. : 0;
        }
        
        int n() const { return 2*l_ + 1;};
        
        std::map<std::string, io::rvec> const& qns() const { return qns_;};
    private:
        int const l_;
        std::map<std::string, io::rvec> qns_;
    };
    
    
    struct SphericalReal {
        SphericalReal() = delete;
        SphericalReal(int l) : l_(l) {
            io::rvec temp(4*l_ + 2, 1.);
            
            qns_["N"] = temp;
            
            for(int i =       0; i < 2*l + 1; ++i) temp[i] =  .5;
            for(int i = 2*l + 1; i < 4*l + 2; ++i) temp[i] = -.5;
            qns_["Sz"] = temp;
        }
        
        std::complex<double> operator()(int os, int m, int s) const {
            if(os/(2*l_ + 1) != s) return .0;
            
            int o = os%(2*l_ + 1);
            
            o -= l_;
            m -= l_;
            
            if(o == 0 && m == 0) return 1.;
            
            std::complex<double> const I(.0, 1.);
            
            if(o < 0) {
                if(m ==  o) return                 I/std::sqrt(2.);
                if(m == -o) return -(o%2 ? -I  : I )/std::sqrt(2.);
            } else {
                if(m == -o) return                1./std::sqrt(2.);
                if(m ==  o) return  (o%2 ? -1. : 1.)/std::sqrt(2.);
            }
            
            return .0;
        };
        
        int n() const { return 2*l_ + 1;};
        
        std::map<std::string, io::rvec> const& qns() const { return qns_;};
    private:
        int const l_;
        std::map<std::string, io::rvec> qns_;
    };
    
    
    struct SphericalCoupled {
        SphericalCoupled() = delete;
        SphericalCoupled(int l) : l_(l) {
            io::rvec temp(4*l_ + 2, 1.);
            
            qns_["N"] = temp;

            for(int i =    0; i <      2*l_; ++i) temp[i] = i -   l_ + .5;
            for(int i = 2*l_; i <  4*l_ + 2; ++i) temp[i] = i - 3*l_ - .5;
            qns_["Jz"] = temp;
        }
        ~SphericalCoupled() = default;
        
        std::complex<double> operator()(int jm_j, int m, int s) const {
            double j   = jm_j < 2*l_ ?        l_ - .5 :          l_ + .5;
            double m_j = jm_j < 2*l_ ? jm_j - l_ + .5 : jm_j - 3*l_ - .5;
            
            return cg::clebschGordan(l_, m - l_, 0.5, s - 0.5, j, m_j); // m s oder s m ????
        };
        
        int n() const { return 2*l_ + 1;};
        
        std::map<std::string, io::rvec> const& qns() const { return qns_;};
    private:
        int const l_;
        std::map<std::string, io::rvec> qns_;
    };
    
    
    struct Transformation {
        Transformation(int const N, jsx::value const& jTransformation) :
        I_(jTransformation.is<jsx::empty_t>() ? N : jTransformation.size()),
        J_(N),
        t_(I_*J_, .0){
            if(jTransformation.is<jsx::empty_t>())
                for(int f = 0; f < N; ++f) (*this)(f, f) = 1.;
            else {
                for(int f = 0; f < I_; ++f)
                    if(jTransformation(f).size() != J_)
                        throw std::runtime_error("opt: invalid transformation matrix");
                
                for(int f_new = 0; f_new < I_; ++f_new)
                    for(int f_old = 0; f_old < J_; ++f_old)
                        (*this)(f_new, f_old) = jTransformation(f_new)(f_old).real64();
            }
            
        }
        
        double operator()(int i, int j) const {
            return t_[i*J_ + j];
        };
        
        int I() const { return I_;};
        int J() const { return J_;};
        
    private:
        int const I_, J_;
        std::vector<double> t_;
        
        double& operator()(int i, int j) {
            return t_[i*J_ + j];
        };
    };
    
    struct Basis {
        template<typename Type>
        Basis(Type const& type, jsx::value const& jTransformation) :
        n_(type.n()),
        transformation_(2*n_, jTransformation) {

            u_.resize(transformation_.I()*transformation_.J(), .0);
            for(int f_new = 0; f_new < transformation_.I(); ++f_new)
                for(int f_old = 0; f_old < transformation_.J(); ++f_old)
                    for(int m = 0; m < n_; ++m)
                        for(int s = 0; s < 2; ++s)
                            (*this)(f_new, m, s) += transformation_(f_new, f_old)*type(f_old, m, s);
            
             for(auto const& qn : type.qns()) {
                 auto const& qn_old = qn.second; io::rvec qn_new;
                 
                 for(int f_new = 0; f_new < transformation_.I(); ++f_new) {
                     std::set<double> count;
                     
                     for(int f_old = 0; f_old < transformation_.J(); ++f_old)
                         if(transformation_(f_new, f_old) != .0) count.insert(qn_old.at(f_old));
                     
                     if(count.size() == 1) qn_new.push_back(*count.begin());
                 }
                 
                 if(qn_new.size() == static_cast<std::size_t>(transformation_.I())) qns_[qn.first] = qn_new;
             }
        }

        int N() const { return transformation_.I();};
        int n() const { return n_;};
    
        std::complex<double> operator()(int f, int m, int s) const {
            return u_[f*2*n_ + 2*m + s];
        };
        
        std::map<std::string, io::rvec> const& qns() const { return qns_;};
    private:
        int const n_;
        Transformation const transformation_;
        std::vector<std::complex<double>> u_;
        std::map<std::string, io::rvec> qns_;
        
        std::complex<double>& operator()(int f, int m, int s) {
            return u_[f*2*n_ + 2*m + s];
        };
    };
    
    
    struct Generic {
        Generic() = delete;
        Generic(int n) : n_(n) {
            io::rvec temp(2*n_, 1.);
            
            qns_["N"] = temp;
            
            for(int i =  0; i   < n_; ++i) temp[i] =  .5;
            for(int i = n_; i < 2*n_; ++i) temp[i] = -.5;
            qns_["Sz"] = temp;
        };
        ~Generic() = default;
        
        std::complex<double> operator()(int ms, int m, int s) const {
            return ( ms%n_ == m )&&( ms/n_ == s ) ? 1. : 0;
        }
        
        int n() const { return n_;};
        int N() const { return 2*n_;};
        
        std::map<std::string, io::rvec> const& qns() const { return qns_;};
    private:
        int const n_;
        std::map<std::string, io::rvec> qns_;
    };
};

#endif






