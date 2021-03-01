#ifndef INCLUDE_OPTIONS_OPTIONS_H
#define INCLUDE_OPTIONS_OPTIONS_H

#include "Basis.h"
#include "Interaction.h"
#include "Observable.h"

#include "../JsonX.h"
#include "../io/Vector.h"
#include "../../ctqmc/include/Utilities.h"

#include <ctime>


namespace opt {

    template<typename Value>
    jsx::value transform(Basis<Value> const& basis, Interaction const& interaction) {
        if(interaction.n() != basis.n())
            throw std::runtime_error("opt::get_interaction_tensor: missmatch between interaction and basis dimension");
        
        io::Vector<Value> two_body(basis.N()*basis.N()*basis.N()*basis.N());
        
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
                                                entry += ut::to_val<Value>(ut::conj(basis(f1Dagg, m1, s)*basis(f2Dagg, m2, sp))*
                                                                           interaction(m1, m2, m3, m4)*
                                                                           basis(f1, m3, sp)*basis(f2, m4, s));
                    }
        
        return two_body;
    };
    
    
    // hack !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    template<typename Value>
    inline void complete_hloc(jsx::value& jParams) {
        if(jParams("hloc")("two body").is<jsx::object_t>() && !(jParams("hloc")("two body").is("real") && jParams("hloc")("two body").is("imag") && jParams("hloc")("two body").size() == 2))
            jParams("hloc")("two body") = transform(get_basis<Value>(jParams("basis")), get_interaction(jParams("basis"), jParams("hloc")("two body")));
    };
    
    
    template<typename Value>
    inline void complete_qn(jsx::value const& jParams, jsx::value& jqn) {
        if(jqn.is<jsx::empty_t>()) jqn = jsx::object_t();
        if(!jqn.is("N")) jqn["N"] = jsx::array_t(jParams("hybridisation")("matrix").size(), jsx::real64_t{1.});
        
        std::vector<std::string> qn_to_delete;
        for(auto& qn : jqn.object()){
            if(qn.second.is<jsx::object_t>()) {
                Basis<Value> basis = get_basis<Value>(jParams("basis"));
                
                if(!basis.qns().count(qn.first))
                    throw std::runtime_error("opt: quantum number " + qn.first + " not defined"); // more sophisticated error message !
                
                qn.second = basis.qns().at(qn.first);
            }
            
            auto Qn = ga::construct_sector_qn(jParams("hloc"), qn.second, false);
            if (!Qn.size()){
                mpi::cout << " ... " << qn.first << " does not survive ... " << std::flush;
                qn_to_delete.push_back(qn.first);
            }
        }
        
        for (auto const& key : qn_to_delete){
            auto& qn = jqn.object();
            auto const it = qn.find(key);
            qn.erase(it);
        }
        
    };
    
    
    template<typename Value>
    jsx::value transform(Observable const& tensor, Transformation<Value> const& transformation, bool ising) {
        if(transformation.J() != tensor.N())
            throw std::runtime_error("opt::get_tensor: missmatch between tensor and transformation dimension");
        
        jsx::value jTensors;
        
        io::PrettyMatrix<Value> one_body(transformation.I(), transformation.I());
        for(int fDagg = 0; fDagg < transformation.I(); ++fDagg)
            for(int f = 0; f < transformation.I(); ++f) {
                Value temp = .0;
                
                for(int gDagg = 0; gDagg < transformation.J(); ++gDagg)
                    for(int g = 0; g < transformation.J(); ++g)
                        temp += ut::conj(transformation(fDagg, gDagg))*tensor.t(gDagg, g)*transformation(f, g);
                
                one_body(fDagg, f) = temp;  
            }
        
        jTensors["one body"] = std::move(one_body);
        
        io::Vector<Value> two_body(transformation.I()*transformation.I()*transformation.I()*transformation.I());
        for(int f1Dagg = 0; f1Dagg < transformation.I(); ++f1Dagg)
            for(int f1 = 0; f1 < transformation.I(); ++f1)
                for(int f2Dagg = 0; f2Dagg < transformation.I(); ++f2Dagg)
                    for(int f2 = 0; f2 < transformation.I(); ++f2) {
                        Value temp = .0;
                        
                        if(ising)
                            if(!((f1Dagg == f2 && f2Dagg == f1) || (f1Dagg == f1 && f2Dagg == f2))) continue;
                        
                        for(int g1Dagg = 0; g1Dagg < transformation.J(); ++g1Dagg)
                            for(int g1 = 0; g1 < transformation.J(); ++g1)
                                for(int g2Dagg = 0; g2Dagg < transformation.J(); ++g2Dagg)
                                    for(int g2 = 0; g2 < transformation.J(); ++g2)
                                        temp += ut::conj(transformation(f1Dagg, g1Dagg))*transformation(f1, g1)*
                                                tensor.V(g1Dagg, g1, g2Dagg, g2)*
                                                ut::conj(transformation(f2Dagg, g2Dagg))*transformation(f2, g2);
                        
                        two_body[transformation.I()*transformation.I()*transformation.I()*f1Dagg +
                                 transformation.I()*transformation.I()*f1 +
                                 transformation.I()*f2Dagg +
                                 f2] = temp;
                    }
        
        jTensors["two body"] = two_body;
        
        return jTensors;
    }

    template <typename Value>
    jsx::value read_observable(jsx::value const& jObservable, int const N, bool ising){
        
        jsx::value jTensors;
        
        io::PrettyMatrix<Value> one_body(N,N);
        io::Vector<Value> two_body(N*N*N*N,0);
        
        
        if (jObservable.is("one body")){

            mpi::cout << "one body part ... " << std::flush;
            
            auto tmp = jObservable("one body");
            auto const& one_body_in = jsx::at<io::PrettyMatrix<Value>>(tmp);
            
            if (one_body_in.I() != one_body.I() || one_body_in.J() != one_body.J())
                throw std::runtime_error("one-body matrix is the wrong size");
            
            one_body = std::move(one_body_in);
            
            mpi::cout << " Ok ... " << std::flush;
                
        } else {
            
            mpi::cout << "no one body part found -- setting to 0 ...";
            
        }
        
        jTensors["one body"] = std::move(one_body);
        
        
        
        
        
        if (jObservable.is("two body")){
            
            mpi::cout << "two body part ... " << std::flush;
            
            auto const& two_body_in = jsx::at<io::Vector<Value>>(jObservable("one body"));
        
            if (two_body_in.size() != two_body.size())
                throw std::runtime_error("two-body matrix is the wrong size");
             
            for(int f1Dagg = 0; f1Dagg < N; ++f1Dagg)
                for(int f1 = 0; f1 < N; ++f1)
                    for(int f2Dagg = 0; f2Dagg < N; ++f2Dagg)
                        for(int f2 = 0; f2 < N; ++f2){
                            if (ising)
                                if(!((f1Dagg == f2 && f2Dagg == f1) || (f1Dagg == f1 && f2Dagg == f2))) continue;
                                
                            std::size_t index = N*N*N*f1Dagg + N*N*f1 +  N*f2Dagg + f2;
                                
                            two_body[index] = two_body_in[index];
                        }
        } else {
                   
            mpi::cout << "no two body part found -- setting to 0 ...";
                   
        }
        
        jTensors["two body"] = two_body;
        
        return jTensors;
        
    }
    
    
    template<typename Value>
    inline void complete_observables(jsx::value const& jParams, jsx::value& jObservables, bool ising) {
        
        if(jObservables.is<jsx::empty_t>()) jObservables = jsx::object_t();
        
        std::vector<std::string> obs_to_delete;
        for(auto& obs : jObservables.object()){
            if(!obs.second.size()) {
                Observable observable = get_observable(jParams("basis"), obs.first);
                Transformation<Value> transformation(observable.N(), jParams("basis").is("transformation") ? jParams("basis")("transformation") : jsx::empty_t());
                
                obs.second = transform(observable, transformation, ising);
            } else {
                mpi::cout << "Reading in and checking observable " << obs.first << " ... " << std::flush;
                obs.second = read_observable<Value>(obs.second, jParams("hybridisation")("matrix").size(), ising);
            }

            obs.second = ga::construct_observable<Value>(jParams("hloc"), obs.second, false);
            
            if (obs.second.is<jsx::empty_t>())
                obs_to_delete.push_back(obs.first);
            
            mpi::cout << " Ok" << std::endl;
        }
        
        //clean up observables which dont survive, i.e., meet commutation requirements
        for (auto const& key : obs_to_delete){
            auto& obs = jObservables.object();
            auto const it = obs.find(key);
            obs.erase(it);
        }
        
    }
    
};

#endif //OPTIONS


