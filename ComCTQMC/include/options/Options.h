#ifndef OPT_OPTIONS_H
#define OPT_OPTIONS_H

#include "Basis.h"
#include "Interaction.h"
#include "Observables.h"

#include "../JsonX.h"
#include "../io/Vector.h"

#include <ctime>


namespace opt {
    
    //--------------------------------------------------- interaction -----------------------------------------------------------
    
    template<typename Interaction, typename Basis>
    jsx::value get_interaction(Interaction const& interaction, Basis const& basis) {
        if(interaction.n() != basis.n())
            throw std::runtime_error("opt::get_interaction: missmatch between interaction and basis dimension");
        
        io::rvec two_body(basis.N()*basis.N()*basis.N()*basis.N());
        
        for(int f1Dagg = 0; f1Dagg < basis.N(); ++f1Dagg)
            for(int f2Dagg = 0; f2Dagg < basis.N(); ++f2Dagg)
                for(int f1 = 0; f1 < basis.N(); ++f1)
                    for(int f2 = 0; f2 < basis.N(); ++f2) {
                        auto& entry = two_body[basis.N()*basis.N()*basis.N()*f1Dagg +
                                               basis.N()*basis.N()*f2Dagg +
                                               basis.N()*f1 +
                                               f2];
                        entry = .0;
                        
                        if(interaction.approximation() == "ising")
                            if(!((f1Dagg == f2 && f2Dagg == f1) || (f1Dagg == f1 && f2Dagg == f2))) continue;
                       
                        
                        for(int m1 = 0; m1 < interaction.n(); ++m1)
                            for(int m2 = 0; m2 < interaction.n(); ++m2)
                                for(int m3 = 0; m3 < interaction.n(); ++m3)
                                    for(int m4 = 0; m4 < interaction.n(); ++m4)
                                        for(int s = 0; s < 2; ++s)
                                            for(int sp = 0; sp < 2; ++sp)
                                                entry += (
                                                          std::conj(basis(f1Dagg, m1, s)*basis(f2Dagg, m2, sp))*
                                                          interaction(m1, m2, m3, m4)*
                                                          basis(f1, m3, sp)*basis(f2, m4, s)).real();
                    }
        
        return std::move(two_body);
    };
    
    
    inline void complete_hloc(jsx::value& jParams) {
        if(jParams("hloc")("two body").is<jsx::object_t>()) {
            if(jParams("basis")("orbitals").is<jsx::string_t>()) {
                std::map<std::string, int> lOption{{"s", 0}, {"p", 1}, {"d", 2}, {"f", 3}};
                
                if(!lOption.count(jParams("basis")("orbitals").string()))
                    throw std::runtime_error("opt: orbitals " + jParams("basis")("orbitals").string() + " not defined");
                
                int const l = lOption.at(jParams("basis")("orbitals").string());
                
                std::unique_ptr<Basis> basis;
                
                jsx::value jTransformation = jParams("basis").is("transformation") ? jParams("basis")("transformation") : jsx::empty_t();
                if(jParams("basis")("type").string() == "product real" || jParams("basis")("type").string() == "product")
                    basis.reset(new Basis(SphericalReal(l), jTransformation));
                else if(jParams("basis")("type").string() == "product imag")
                    basis.reset(new Basis(SphericalImag(l), jTransformation));
                else if(jParams("basis")("type").string() == "coupled")
                    basis.reset(new Basis(SphericalCoupled(l), jTransformation));
                else
                    throw std::runtime_error("opt: type " + jParams("basis")("type").string() + " not defined for spherical harmonics basis");
                
                if(jParams("hloc")("two body")("parametrisation").string() == "slater-condon")
                    jParams("hloc")("two body") = get_interaction(SlaterCondon(l, jParams("hloc")("two body")), *basis);
                else
                    throw std::runtime_error("opt: parametrisation " + jParams("hloc")("two body")("parametrisation").string() + " not defined for spherical harmonics");
                
            } else if(jParams("basis")("orbitals").is<jsx::int64_t>()) {
                
                int const n = jParams("basis")("orbitals").int64();
                
                if(jParams("hloc")("two body")("parametrisation").string() == "kanamori")
                    jParams("hloc")("two body") = get_interaction(Kanamori(n, jParams("hloc")("two body")), Generic(n));
                
                if(jParams("hloc")("two body").is<jsx::object_t>())
                    throw std::runtime_error("opt: parametrisation " + jParams("hloc")("two body")("parametrisation").string() + " not defined for generic basis");
                
            } else
                throw std::runtime_error("opt: orbitals not defined");
            
            
        }
    };

    //--------------------------------------------------- quantum numbers -----------------------------------------------------------
    
    inline void complete_qn(jsx::value& jParams) {
        if(!jParams.is("quantum numbers")) jParams["quantum numbers"] = jsx::object_t();
        if(!jParams("quantum numbers").is("N")) jParams("quantum numbers")["N"] = jsx::object_t();
        
        for(auto& qn : jParams("quantum numbers").object())
            if(qn.second.is<jsx::object_t>()) {
                if(jParams("basis")("orbitals").is<jsx::string_t>()) {
                    std::map<std::string, int> lOption{{"s", 0}, {"p", 1}, {"d", 2}, {"f", 3}};
                    
                    if(!lOption.count(jParams("basis")("orbitals").string()))
                        throw std::runtime_error("opt: orbitals " + jParams("basis")("orbitals").string() + " not defined");
                    
                    int const l = lOption.at(jParams("basis")("orbitals").string());
                    
                    std::unique_ptr<Basis> basis;
                    
                    jsx::value jTransformation = jParams("basis").is("transformation") ? jParams("basis")("transformation") : jsx::empty_t();
                    if(jParams("basis")("type").string() == "product real" || jParams("basis")("type").string() == "product")
                        basis.reset(new Basis(SphericalReal(l), jTransformation));
                    else if(jParams("basis")("type").string() == "product imag")
                        basis.reset(new Basis(SphericalImag(l), jTransformation));
                    else if(jParams("basis")("type").string() == "coupled")
                        basis.reset(new Basis(SphericalCoupled(l), jTransformation));
                    else
                        throw std::runtime_error("opt: type " + jParams("basis")("type").string() + " not defined for spherical harmonics basis");
                    
                    if(!basis->qns().count(qn.first))
                        throw std::runtime_error("opt: quantum number " + qn.first + " not defined"); // more sophisticated error message !
                    
                    qn.second = basis->qns().at(qn.first);
                    
                } else if(jParams("basis")("orbitals").is<jsx::int64_t>()) {
                    
                    int const n = jParams("basis")("orbitals").int64();
                    
                    Generic basis(n);
                    
                    if(!basis.qns().count(qn.first))
                        throw std::runtime_error("opt: quantum number " + qn.first + " not defined"); // more sophisticated error message !
                    
                    qn.second = basis.qns().at(qn.first);
                    
                } else
                    throw std::runtime_error("opt: orbitals not defined");
            }
    };
    
    //--------------------------------------------------- observables -----------------------------------------------------------
    
    template<typename Tensor, typename Transformation>
    jsx::value get_observable(Tensor const& tensor, Transformation const& transformation) {
        if(transformation.J() != tensor.N())
            throw std::runtime_error("opt::get_tensor: missmatch between tensor and transformation dimension");
        
        jsx::value jTensors;
        
        jsx::value one_body = jsx::array_t(transformation.I(), jsx::array_t(transformation.I()));
        for(int fDagg = 0; fDagg < transformation.I(); ++fDagg)
            for(int f = 0; f < transformation.I(); ++f) {
                double temp = .0;
                
                for(int gDagg = 0; gDagg < transformation.J(); ++gDagg)
                    for(int g = 0; g < transformation.J(); ++g)
                        temp += transformation(fDagg, gDagg)*tensor.t(gDagg, g)*transformation(f, g);
                
                one_body(fDagg)(f) = temp;
            }
        jTensors["one body"] = std::move(one_body);
        
        io::rvec two_body(transformation.I()*transformation.I()*transformation.I()*transformation.I());
        for(int f1Dagg = 0; f1Dagg < transformation.I(); ++f1Dagg)
            for(int f1 = 0; f1 < transformation.I(); ++f1)
                for(int f2Dagg = 0; f2Dagg < transformation.I(); ++f2Dagg)
                    for(int f2 = 0; f2 < transformation.I(); ++f2) {
                        double temp = .0;
                        
                        for(int g1Dagg = 0; g1Dagg < transformation.J(); ++g1Dagg)
                            for(int g1 = 0; g1 < transformation.J(); ++g1)
                                for(int g2Dagg = 0; g2Dagg < transformation.J(); ++g2Dagg)
                                    for(int g2 = 0; g2 < transformation.J(); ++g2)
                                        temp += transformation(f1Dagg, g1Dagg)*transformation(f1, g1)*
                                                tensor.V(g1Dagg, g1, g2Dagg, g2)*
                                                transformation(f2Dagg, g2Dagg)*transformation(f2, g2);
                        
                        two_body[transformation.I()*transformation.I()*transformation.I()*f1Dagg +
                                 transformation.I()*transformation.I()*f1 +
                                 transformation.I()*f2Dagg +
                                 f2] = temp;
                    }
        jTensors["two body"] = two_body;
        
        return jTensors;
    }
    
    inline void complete_observables(jsx::value& jParams) {
        if(!jParams.is("observables")) jParams["observables"] = jsx::object_t();
        
        for(auto& obs : jParams("observables").object())
            if(!obs.second.size()) {
                if(jParams("basis")("orbitals").is<jsx::string_t>()) {
                    std::map<std::string, int> lOption{{"s", 0}, {"p", 1}, {"d", 2}, {"f", 3}};
                    
                    if(!lOption.count(jParams("basis")("orbitals").string()))
                        throw std::runtime_error("opt: orbitals " + jParams("basis")("orbitals").string() + " not defined");
                    
                    int const l = lOption.at(jParams("basis")("orbitals").string());
                    
                    Transformation transformation(4*l + 2, jParams("basis").is("transformation") ? jParams("basis")("transformation") : jsx::empty_t());

                    if(jParams("basis")("type").string() == "product real" || jParams("basis")("type").string() == "product") {
                        if(obs.first == "S2")
                            obs.second = get_observable(S2(2*l + 1), transformation);
                    } else if(jParams("basis")("type").string() == "product imag") {
                        if(obs.first == "S2")
                            obs.second = get_observable(S2(2*l + 1), transformation);
                        /*else if(obs.first == "L2")
                            obs.second = get_observable(L2(l), transformation);*/
                    } else if(jParams("basis")("type").string() == "coupled") {
                        if(obs.first == "J2")
                            obs.second = get_observable(J2(l), transformation);
                    } else
                        throw std::runtime_error("opt: type " + jParams("basis")("type").string() + " not defined for spherical harmonics basis");
                    
                    if(!obs.second.size())
                        throw std::runtime_error("opt: observable " + obs.first + " not defined for spherical harmonics of " + jParams("basis")("type").string() + " type");
                    
                } else if(jParams("basis")("orbitals").is<jsx::int64_t>()) {
                    
                    int const n = jParams("basis")("orbitals").int64();
                    
                    Transformation transformation(2*n, jsx::empty_t());
                    
                    if(obs.first == "S2")
                        obs.second = get_observable(S2(n), transformation);
                    
                    if(!obs.second.size())
                        throw std::runtime_error("opt: observable " + obs.first + " not defined for generic orbitals");
                    
                } else
                    throw std::runtime_error("opt: orbitals not defined");
            }
    }
};

#endif //OPTIONS


