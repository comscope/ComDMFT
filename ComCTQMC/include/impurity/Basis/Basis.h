#ifndef BASIS
#define BASIS

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

#include "../../JsonX.h"
#include "../../IO.h"

//////////////////////////////////////////////////
//Scheisse ordnig:   m,s oder s,m ???
//////////////////////////////////////////////////

namespace Basis {

    struct Basis {
        
        int dim() const { return dim_;};
        int l() const { return l_;};
        virtual std::string type() const = 0;
        
        std::complex<double> operator()(int f, int m, int s) const {
            return tu_[2*(2*l() + 1)*f + 2*m + s];
        };
        std::map<std::string, io::rvec> const& qns() const {
            return qns_;
        };
        
        virtual ~Basis() = default;
        
    protected:
        std::map<std::string, io::rvec> qns_;
        
        Basis(jsx::value const& jBasis) {
            auto const shell = jBasis("shell").string();
            
            if(shell == "s")
                l_ = 0;
            else if(shell == "p")
                l_ = 1;
            else if(shell == "d")
                l_ = 2;
            else if(shell == "f")
                l_ = 3;
            else
                throw std::runtime_error("shell " + shell + " not defined.");
            
            dim_ = jBasis.is("transformation") ? jBasis("transformation").size() : 2*(2*l_ + 1);
            
            t_.resize(dim_*2*(2*l_ + 1), .0);
            tu_.resize(dim_*2*(2*l_ + 1), .0);
            
            if(jBasis.is("transformation")) {
                for(int f = 0; f < dim_; ++f)
                    if(jBasis("transformation")(f).size() != static_cast<std::size_t>(2*(2*l_ + 1)))
                        throw std::runtime_error("Transformation matrix has wrong size");
                
                for(int f_new = 0; f_new < dim_; ++f_new)
                    for(int f_old = 0; f_old < 2*(2*l_ + 1); ++f_old)
                       T(f_new, f_old) = jBasis("transformation")(f_new)(f_old).real64();
            } else
                for(int f = 0; f < 2*(2*l_ + 1); ++f)
                    T(f, f) = 1.;
            
            //test if transformation is unitary
        };

        double& T(int f_new, int f_old) {
            return t_[2*(2*l_ + 1)*f_new + f_old];
        };
        
        std::complex<double>& TU(int f, int m, int s) {
            return tu_[2*(2*l() + 1)*f + 2*m + s];
        };
        
    private:
        int l_;
        int dim_;
        std::vector<double> t_;
        std::vector<std::complex<double>> tu_;
        
        Basis() = delete;
        Basis(Basis const&) = delete;
        Basis(Basis&&) = delete;
        Basis& operator=(Basis const&) = delete;
        Basis& operator=(Basis&&) = delete;
    };
    
    
    struct ProductBasis : Basis {
        ProductBasis() = delete;
        ProductBasis(jsx::value const& jBasis) :
        Basis(jBasis) {
            for(int f = 0; f < dim(); ++f)
                for(int m = 0; m < 2*l() + 1; ++m)
                    for(int s = 0; s < 2; ++s)
                        for(int o = 0; o < 2*l() + 1; ++o)
                            TU(f, m, s)  += T(f, (2*l() + 1)*s + o)*U(o, m);
            
            qns_["Sz"].resize(dim());
            for(int f = 0; f < dim(); ++f) {
                bool isUp = false;
                for(int o = 0; o < 5; ++o)
                    if(T(f,               o) != .0) isUp = true;

                bool isDown = false;
                for(int o = 0; o < 5; ++o)
                    if(T(f, (2*l() + 1) + o) != .0) isDown = true;
                
                if(qns_.at("Sz").size() && isUp != isDown)
                    qns_.at("Sz").at(f) = isUp ? .5 : -.5;
                else
                    qns_.at("Sz").clear();
            }
        };
        
        std::string type() const {
            return "product";
        };
        
        ~ProductBasis() = default;
    private:
        std::complex<double> U(int o, int m) const {
            o -= l();
            m -= l();
            
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
    };
    
    
    struct CoupledBasis : Basis {
        CoupledBasis() = delete;
        CoupledBasis(jsx::value const& jBasis) :
        Basis(jBasis) {
            for(int f = 0; f < dim(); ++f)
                for(int m = 0; m < 2*l() + 1; ++m)
                    for(int s = 0; s < 2; ++s)
                        for(int jm_j = 0; jm_j < 2*(2*l() + 1); ++jm_j)
                            TU(f, m, s)  += T(f, jm_j)*U(jm_j, m, s);
            
            qns_["Jz"].resize(dim());
            for(int f = 0; f < dim(); ++f) {
                std::vector<double> temp;

                if(    T(f,        2*l()    ) != .0                   ) temp.push_back(     -l() - 0.5);
                
                for(int m_j = 0; m_j < 2*l(); ++m_j)
                    if(T(f,  m_j + 2*l() + 1) != .0 || T(f, m_j) != .0) temp.push_back(m_j - l() + 0.5);
                
                if(    T(f,        4*l() + 1) != .0                   ) temp.push_back(      l() + 0.5);
                    
                if(qns_.at("Jz").size() && temp.size() == 1)
                    qns_.at("Jz").at(f) = temp.at(0);
                else
                    qns_.at("Jz").clear();
            }
        };
        
        std::string type() const {
            return "coupled";
        };
        
        ~CoupledBasis() = default;
    private:
        double U(int jm_j, int m, int s) const {
            double j   = jm_j < 6 ?        2.5 :        3.5;
            double m_j = jm_j < 6 ? jm_j - 2.5 : jm_j - 9.5;
            
            return cg::clebschGordan(0.5, s - 0.5, l(), m - l(), j, m_j); // m s oder s m ????
        };
    };
    
    
    std::unique_ptr<Basis> get_basis(jsx::value const& jBasis) {
        std::unique_ptr<Basis> basis;
        
        std::cout << "Reading basis ... ";
        
        if(jBasis("type").string() == "product")
            basis = std::unique_ptr<Basis>(new ProductBasis(jBasis));
        else if(jBasis("type").string() == "coupled")
            basis = std::unique_ptr<Basis>(new CoupledBasis(jBasis));
        else
            throw std::runtime_error("type " + jBasis("type").string() + " not defined.");
        
        std::cout << "Ok" << std::endl;
        
        return basis;
    };
};

#endif






