#ifndef INCLUDE_OPTIONS_OBSERVABLE_H
#define INCLUDE_OPTIONS_OBSERVABLE_H


#include "sphericalharmonics/Observables.h"
#include "model/Observables.h"

#include "../JsonX.h"
#include "../io/Vector.h"


namespace opt {

    struct Observable {
        Observable() = delete;
        template<typename Tensor>
        Observable(Tensor const& tensor) :
        N_(tensor.N()), one_body_(N_*N_, .0), two_body_(N_*N_*N_*N_, .0) {
            for(int fDagg = 0; fDagg < N_; ++fDagg)
                for(int f = 0; f < N_; ++f)
                    one_body_[N_*fDagg + f] = tensor.t(fDagg, f);
            
            for(int f1Dagg = 0; f1Dagg < N_; ++f1Dagg)
                for(int f1 = 0; f1 < N_; ++f1)
                    for(int f2Dagg = 0; f2Dagg < N_; ++f2Dagg)
                        for(int f2 = 0; f2 < N_; ++f2)
                            two_body_[N_*N_*N_*f1Dagg +
                                      N_*N_*f1 +
                                      N_*f2Dagg +
                                      f2] = tensor.V(f1Dagg, f1, f2Dagg, f2);
        }
        Observable(Observable const&) = delete;
        Observable(Observable&& ) = default;
        Observable& operator=(Observable const& ) = delete;
        Observable& operator=(Observable&& ) = default;
        ~Observable() = default;

        int N() const {
            return N_;
        };
        
        double t(int fDagg, int f) const {
            return one_body_[N_*fDagg + f];
        };
        
        double V(int f1Dagg, int f1, int f2Dagg, int f2) const {
            return two_body_[N_*N_*N_*f1Dagg +
                             N_*N_*f1 +
                             N_*f2Dagg +
                             f2];
        };
        
    private:
        int N_;
        std::vector<double> one_body_;
        std::vector<double> two_body_;
    };
    
    
    Observable get_observable(jsx::value const& jBasis, std::string const name) {
        if(jBasis("orbitals").is<jsx::string_t>()) {

            if(jBasis("type").string() == "product real" || jBasis("type").string() == "product") {
                if(name == "S2")
                    return Observable(sphericalharmonics::S2(sphericalharmonics::Basis(jBasis)));
            } else if(jBasis("type").string() == "product imag") {
                if(name == "S2")
                    return Observable(sphericalharmonics::S2(sphericalharmonics::Basis(jBasis)));
                //else if(obs.first == "L2")
                //    return Observable(sphericalharmonics::L2(l), transformation);
            } else if(jBasis("type").string() == "coupled") {
                if(name == "J2")
                    return Observable(sphericalharmonics::J2(sphericalharmonics::Basis(jBasis)));
            }
            
            throw std::runtime_error("opt: observable " + name + " not defined for spherical harmonics");

        } else if(jBasis("orbitals").is<jsx::int64_t>()) {

            if(name == "S2")
                return Observable(model::S2(model::Basis(jBasis)));
            
            throw std::runtime_error("opt: observable " + name + " not defined for model basis");
            
        } else
            
            throw std::runtime_error("opt: orbitals not defined");
    }
    
};

#endif //OPTIONS


