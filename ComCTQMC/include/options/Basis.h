#ifndef INCLUDE_OPTIONS_BASIS_H
#define INCLUDE_OPTIONS_BASIS_H

#include "Transformation.h"
#include "sphericalharmonics/Basis.h"
#include "model/Basis.h"

#include "../JsonX.h"
#include "../io/Vector.h"


namespace opt {
    
    template<typename Value>
    struct Basis {
        Basis() = delete;
        template<typename Type>
        Basis(Type const& type, jsx::value const& jTransformation) :
        n_(type.n()), transformation_(2*n_, jTransformation) {
            
            u_.resize(transformation_.I()*transformation_.J(), .0);
            for(int f_new = 0; f_new < transformation_.I(); ++f_new)
                for(int f_old = 0; f_old < transformation_.J(); ++f_old)
                    for(int m = 0; m < n_; ++m)
                        for(int s = 0; s < 2; ++s)
                            u_[f_new*2*n_ + 2*m + s] += transformation_(f_new, f_old)*type(f_old, m, s);

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
        Basis(Basis const&) = delete;
        Basis(Basis&& ) = default;
        Basis& operator=(Basis const& ) = delete;
        Basis& operator=(Basis&& ) = default;
        ~Basis() = default;
        
        int n() const {
            return n_;
        };
        
        int N() const {
            return transformation_.I();
        };
        
        std::complex<double> operator()(int f, int m, int s) const {
            return u_[f*2*n_ + 2*m + s];
        };
        
        std::map<std::string, io::rvec> const& qns() const {
            return qns_;
        };
        
    private:
        int n_;
        Transformation<Value> const transformation_;
        std::vector<std::complex<double>> u_;
        std::map<std::string, io::rvec> qns_;
    };
    
    
    template<typename Value>
    Basis<Value> get_basis(jsx::value jBasis) {
        jsx::value jTransformation = jBasis.is("transformation") ? jBasis("transformation") : jsx::empty_t();

        if(jBasis("orbitals").is<jsx::string_t>()) {

            if(jBasis("type").string() == "product real" || jBasis("type").string() == "product")
                return Basis<Value>(sphericalharmonics::Real(jBasis), jTransformation);
            else if(jBasis("type").string() == "product imag")
                return Basis<Value>(sphericalharmonics::Imag(jBasis), jTransformation);
            else if(jBasis("type").string() == "coupled")
                return Basis<Value>(sphericalharmonics::Coupled(jBasis), jTransformation);
            else
                throw std::runtime_error("opt: type " + jBasis("type").string() + " not defined for spherical harmonics basis");

        } else if(jBasis("orbitals").is<jsx::int64_t>()) {
            
            return Basis<Value>(model::Basis(jBasis), jTransformation);
            
        } else
            
            throw std::runtime_error("opt: orbitals not defined");
    }
    
};

#endif //OPTIONS


