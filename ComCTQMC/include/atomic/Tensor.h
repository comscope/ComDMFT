#ifndef INCLUDE_ATOMIC_TENSOR_H
#define INCLUDE_ATOMIC_TENSOR_H

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
#include "../io/Matrix.h"

namespace ga {
    
    template<typename Value>
    struct Tensor {
        struct Interaction {};
        
        Tensor() = delete;
        Tensor(jsx::value jTensor) :
        one_body_(std::move(jsx::at<io::PrettyMatrix<Value>>(jTensor("one body")))),
        two_body_(std::move(jsx::at<io::Vector<Value>>(jTensor("two body")))),
        N_(one_body_.I()) {
            if(one_body_.I() != one_body_.J())
                throw std::runtime_error("ga::Tensor: one body matrix is not square");
            
            if(N_*N_*N_*N_ != static_cast<int>(two_body_.size()))
                throw std::runtime_error("Tensor: one and two body tensor dimensions not compatible");
            
            for(auto& entry : two_body_) entry /= 2.;
        };
        Tensor(Tensor const& tensor, Interaction) :
        one_body_(tensor.N_, tensor.N_),
        two_body_(tensor.two_body_),
        N_(tensor.N_) {
        };
        Tensor(Tensor const&) = delete;
        Tensor(Tensor&&) = default;
        Tensor& operator=(Tensor const& rhs) = delete;
        Tensor& operator=(Tensor&&) = default;
        
        int N() const { return N_;};
        
        Value t(int f1, int f2) const {
            return one_body_(f1, f2);
        };
        
        Value V(int f1, int f2, int f3, int f4) const {
            return two_body_[N_*N_*N_*f1 +
                             N_*N_*f2 +
                             N_*f3 +
                             f4];
        };

    private:
        io::PrettyMatrix<Value> one_body_;
        io::Vector<Value> two_body_;
        
        int const N_;
    };

};

#endif






