#ifndef CTQMC_INCLUDE_IMPURITY_TENSOR_H
#define CTQMC_INCLUDE_IMPURITY_TENSOR_H

#include <cmath>
#include <iostream>
#include <vector>

#include "../Utilities.h"
#include "../../../include/JsonX.h"
#include "../../../include/io/Vector.h"


namespace imp {
    
    struct TensorIndices {
        TensorIndices() = delete;
        TensorIndices(int flavor1, int flavor2, int flavor3) : flavor1_(flavor1), flavor2_(flavor2), flavor3_(flavor3) {};
        
        int flavor1() const { return flavor1_;};
        int flavor2() const { return flavor2_;};
        int flavor3() const { return flavor3_;};
        
    private:
        int flavor1_, flavor2_, flavor3_;
    };
    
    
    // all this needs to be modified in the presence of superconductivity
    
    
    template<typename Value>
    struct Tensor {
        Tensor() = delete;
        Tensor(jsx::value jTwoBody, int N) :
        N_(N),
        tensor_(jsx::at<io::Vector<Value>>(jTwoBody)),
        non_zero_(N_) {
            if(N_*N_*N_*N_ != tensor_.size())
                throw std::runtime_error("imp::Tensor: interaction tensor has wrong size");
            
            io::Vector<Value> temp(tensor_.size(), .0);
            for(int i = 0; i < N_; ++i)
                for(int j = 0; j < N_; ++j)
                    for(int k = 0; k < N_; ++k)
                        for(int l = 0; l < N_; ++l) {
                            auto const& element = ((*this)(i, j, k, l) - (*this)(j, i, k, l) - (*this)(i, j, l, k) + (*this)(j, i, l, k))/4.;
                            
                            if(std::abs(element) > 1.e-14) {
                                temp[N_*N_*N_*i + N_*N_*j + N_*k + l] = element;
                                non_zero_[i].push_back({2*j + 1, 2*k, 2*l});
                            }  
                        }
            
            tensor_ = temp;
        };
        
        Value operator()(int i, int j, int k, int l) const {
            return tensor_[N_*N_*N_*i + N_*N_*j + N_*k + l];
        };
        
        std::vector<TensorIndices> const& nonzero(int i) const{
            return non_zero_[i];
        };
        
        //std::vector<TensorIndices> const& nonzerodagg(int i) const{
        //    return non_zero_[i];
        //};
        
    private:
        int const N_;
        io::Vector<Value> tensor_;
        std::vector<std::vector<TensorIndices>> non_zero_;
        
        Value& operator()(int i, int j, int k, int l) {
            return tensor_[N_*N_*N_*i + N_*N_*j + N_*k + l];
        };
        
    };
    
}


#endif
