#ifndef TENSOR
#define TENSOR

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
#include "../io/Vector.h"

namespace ga {
    
    struct Tensor {
        Tensor() = delete;
        Tensor(jsx::value jTensor) :
        N_(jTensor("one body").size()),
        one_body_(N_*N_, .0),
        two_body_(std::move(jsx::at<io::rvec>(jTensor("two body")))) {
            for(int fDagg = 0; fDagg < N_; ++fDagg) {
                if(jTensor("one body")(fDagg).size() != static_cast<std::size_t>(N_))
                    throw std::runtime_error("Tensor: one body matrix has wrong size");
                
                for(int f = 0; f < N_; ++f)
                    one_body_[N_*fDagg + f] = jTensor("one body")(fDagg)(f).real64();
            }
            
            if(N_*N_*N_*N_ != static_cast<int>(two_body_.size()))
                throw std::runtime_error("Tensor: one and two body tensor dimensions not compatible");
            
            for(auto& entry : two_body_) entry /= 2.;
        };        
        Tensor(Tensor const&) = delete;
        Tensor(Tensor&&) = default;
        Tensor& operator=(Tensor const& rhs) = delete;
        Tensor& operator=(Tensor&&) = default;
        
        int N() const { return N_;};
        
        double t(int f1, int f2) const {
            return one_body_[N_*f1 + f2];
        };
        
        double V(int f1, int f2, int f3, int f4) const {
            return two_body_[N_*N_*N_*f1 +
                             N_*N_*f2 +
                             N_*f3 +
                             f4];
        };

    private:
        int const N_;
        std::vector<double> one_body_;
        io::rvec two_body_;
    };
    
    //-----------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------
    
    struct TwoBody {
        TwoBody(Tensor const& tensor) : tensor_(tensor) {};
        double t(int f1, int f2) const {
            return .0;
        };
        double V(int f1, int f2, int f3, int f4) const {
            return tensor_.V(f1, f2, f3, f4);
        };
        int N() const {
            return tensor_.N();
        };
    private:
        Tensor const& tensor_;
    };
};

#endif






