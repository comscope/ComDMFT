#ifndef INCLUDE_OPTIONS_INTERACTION_H
#define INCLUDE_OPTIONS_INTERACTION_H

#include "sphericalharmonics/Basis.h"
#include "sphericalharmonics/SlaterCondon.h"

#include "model/Basis.h"
#include "model/Kanamori.h"

#include "../JsonX.h"
#include "../io/Vector.h"

namespace opt {
    
    struct Interaction {
        Interaction() = delete;
        template<typename Type>
        Interaction(Type const& type) :
        n_(type.n()),
        approximation_(type.approximation()),
        tensor_(n_*n_*n_*n_) {
            for(int m1 = 0; m1 < n_; ++m1)
                for(int m2 = 0; m2 < n_; ++m2)
                    for(int m3 = 0; m3 < n_; ++m3)
                        for(int m4 = 0; m4 < n_; ++m4)
                            tensor_[n_*n_*n_*m1 +
                                    n_*n_*m2 +
                                    n_*m3 +
                                    m4] = type(m1, m2, m3, m4);
        }
        Interaction(Interaction const& ) = delete;
        Interaction(Interaction&& ) = default;
        Interaction& operator=(Interaction const& ) = delete;
        Interaction& operator=(Interaction&& ) = default;
        ~Interaction() = default;
        
        double operator()(int m1, int m2, int m3, int m4) const {
            return tensor_[n_*n_*n_*m1 +
                           n_*n_*m2 +
                           n_*m3 +
                           m4];
        };
        
        int n() const {
            return n_;
        };
        
        std::string approximation() const {
            return approximation_;
        };
        
    private:
        int const n_;
        std::string const approximation_;
        std::vector<double> tensor_;
    };
    
    
    
    Interaction get_interaction(jsx::value jBasis, jsx::value jTwoBody) {
        if(jBasis("orbitals").is<jsx::string_t>()) {
            
            sphericalharmonics::Basis basis(jBasis);
            
            if(jTwoBody("parametrisation").string() == "slater-condon")
                return Interaction(sphericalharmonics::SlaterCondon(basis, jTwoBody));
            else
                throw std::runtime_error("opt: parametrisation" + jTwoBody("parametrisation").string() + " not defined for spherical harmonics basis");
            
        } else if(jBasis("orbitals").is<jsx::int64_t>()) {
            
            model::Basis basis(jBasis);
            
            if(jTwoBody("parametrisation").string() == "kanamori")
                return Interaction(model::Kanamori(basis, jTwoBody));
            else
                throw std::runtime_error("opt: parametrisation" + jTwoBody("parametrisation").string() + " not defined for model basis");
            
        } else
            throw std::runtime_error("opt: orbitals not defined");
    }
    
};

#endif //OPTIONS


