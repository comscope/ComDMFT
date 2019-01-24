#ifndef OPTIONS
#define OPTIONS

#include "Basis/Basis.h"
#include "Interaction/Interaction.h"

#include "../JsonX.h"
#include "../IO.h"


namespace options {
    
    inline jsx::value get_two_body(Interaction::SlaterCondon const& interaction, Basis::Basis const& basis) {
        int const dim = basis.dim();
        
        io::rvec two_body; two_body.resize(dim*dim*dim*dim, .0);
        for(int f1Dagg = 0; f1Dagg < dim; ++f1Dagg)
            for(int f2Dagg = 0; f2Dagg < dim; ++f2Dagg)
                for(int f1 = 0; f1 < dim; ++f1)
                    for(int f2 = 0; f2 < dim; ++f2) {
                        two_body[dim*dim*dim*f1Dagg +
                                 dim*dim*f2Dagg +
                                 dim*f1 +
                                 f2] = .0;
                        
                        // this needs to be generalized !!!!!!!!!
                        
                        /*if(interaction.approximation() == "truncated") {       
                             if(!((f1Dagg%(dim/2) == f2%(dim/2) && f2Dagg%(dim/2) == f1%(dim/2)) || (f1Dagg%(dim/2) == f2Dagg%(dim/2) && f1%(dim/2) == f2%(dim/2)) || (f1Dagg%(dim/2) == f1%(dim/2) && f2Dagg%(dim/2) == f2%(dim/2))))
                                continue;
                        } else */
                        if(interaction.approximation() == "ising") {
                            if(!((f1Dagg == f2 && f2Dagg == f1) || (f1Dagg == f1 && f2Dagg == f2)))
                                continue;
                        } else if(interaction.approximation() != "none")
                            throw std::runtime_error("approximation " + interaction.approximation() + " not defined.");
                        
                        std::complex<double> temp = .0;
                        
                        for(int m1 = 0; m1 < (2*basis.l() + 1); ++m1)
                            for(int m2 = 0; m2 < (2*basis.l() + 1); ++m2)
                                for(int m3 = 0; m3 < (2*basis.l() + 1); ++m3)
                                    for(int m4 = 0; m4 < (2*basis.l() + 1); ++m4)
                                        for(int s = 0; s < 2; ++s)
                                            for(int sp = 0; sp < 2; ++sp)
                                                temp +=
                                                std::conj(basis(f1Dagg, m1, s)*basis(f2Dagg, m2, sp))*
                                                interaction(m1, m2, m3, m4)*
                                                basis(f1, m3, sp)*basis(f2, m4, s);

                        two_body[dim*dim*dim*f1Dagg +
                                 dim*dim*f2Dagg +
                                 dim*f1 +
                                 f2] = temp.real()/2.;
                    }
        
        return two_body;
    };
    
    
    inline jsx::value get_hloc(jsx::value const& jParams) {
        jsx::value jHloc;
        
        jHloc["one body"] = jParams("hloc")("one body");

        std::unique_ptr<Basis::Basis> basis;
        
        if(jParams("hloc")("two body").type() == jsx::type::object) {
            if(basis.get() == nullptr) basis = Basis::get_basis(jParams("basis"));
            
            if(jParams("hloc")("two body")("parametrisation").string() == "slater-condon") {
                jHloc["two body"] = get_two_body(Interaction::SlaterCondon(basis->l(), jParams("hloc")("two body")), *basis);
            } else // Kanamori
                throw std::runtime_error("parametrisation " + jParams("hloc")("two body")("parametrisation").string() + " not defined.");
        } else
            jHloc["two body"] = jParams("hloc")("two body");
        
        
        if(jParams("hloc").is("quantum numbers"))
            for(auto const& qn : jParams("hloc")("quantum numbers").object()) {
                if(qn.second.type() == jsx::type::object) {
                    if(basis.get() == nullptr) basis = Basis::get_basis(jParams("basis"));
                    
                    auto const& qns = basis->qns();
                    if(qns.count(qn.first))
                        jHloc["quantum numbers"][qn.first] = qns.at(qn.first);
                    else
                        throw std::runtime_error("quantum number " + qn.first + " not defined for basis type " + basis->type());
                } else
                    jHloc["quantum numbers"][qn.first] = qn.second;
            }
        
        return jHloc;
    };
};

#endif //OPTIONS


