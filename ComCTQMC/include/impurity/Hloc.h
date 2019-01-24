#ifndef HLOC
#define HLOC

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

#include "../JsonX.h"
#include "../IO.h"

namespace Ga {
    
    struct Hloc {
        Hloc() = delete;
        
        Hloc(jsx::value jHloc) :
        flavor_number_(jHloc("one body").size()),
        one_body_(flavor_number_*flavor_number_, .0),
        two_body_(std::move(jsx::at<io::rvec>(jHloc("two body")))) {
            
            for(int fDagg = 0; fDagg < flavor_number_; ++fDagg) {
                if(jHloc("one body")(fDagg).size() != static_cast<std::size_t>(flavor_number_))
                    throw std::runtime_error("Hloc: one body matrix has wrong size");
                
                for(int f = 0; f < flavor_number_; ++f)
                    one_body_[flavor_number_*fDagg + f] = jHloc("one body")(fDagg)(f).real64();
            }
            
            if(flavor_number_*flavor_number_*flavor_number_*flavor_number_ != static_cast<int>(two_body_.size()))
                throw std::runtime_error("Hloc: one and two body tensor dimensions not compatible");
            
            if(jHloc.is("quantum numbers"))
                for(auto& qn : jHloc("quantum numbers").object())
                    qns_[qn.first] = jsx::at<io::rvec>(qn.second);
            
            qns_["N"] = std::vector<double>(flavor_number_, 1.);
        };
        
        Hloc(Hloc const&) = delete;
        Hloc(Hloc&&) = delete;
        Hloc& operator=(Hloc const& rhs) = delete;
        Hloc& operator=(Hloc&&) = delete;
        
        int N() const { return flavor_number_;};
        
        double t(int fDagg, int f) const {
            return one_body_[flavor_number_*fDagg + f];
        };
        
        double V(int f1Dagg, int f2Dagg, int f1, int f2) const {
            return two_body_[flavor_number_*flavor_number_*flavor_number_*f1Dagg +
                             flavor_number_*flavor_number_*f2Dagg +
                             flavor_number_*f1 +
                             f2];
        };
        
        std::map<std::string, std::vector<double>> const& qns() const {
            return qns_;
        };
    private:
        int const flavor_number_;
        std::vector<double> one_body_;
        io::rvec two_body_;
        std::map<std::string, std::vector<double>> qns_;
    };
};

#endif






