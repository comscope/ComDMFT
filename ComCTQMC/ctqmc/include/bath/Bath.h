#ifndef BATH_BATH_H
#define BATH_BATH_H

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include "Algebra.h"
#include "Hyb.h"
#include "../Utilities.h"

//Vielleicht insert und erase noch verändern wie in CTBI hürütütütüt pluschter pluschter

namespace bath {
    
    struct Operator {
        Operator(ut::KeyType key, int flavor, double* ptr) : key_(key), flavor_(flavor), ptr_(ptr) {};
        
        ut::KeyType key() const { return key_;};
        int flavor() const { return flavor_;};
        double* ptr() const { return ptr_;};
        
    private:
        ut::KeyType key_;
        int flavor_;
        double* ptr_;
    };
    
    
    namespace itf {
        
        struct Bath {
            virtual std::vector<Operator> const& opsL() const = 0;
            virtual std::vector<Operator> const& opsR() const = 0;
            
            virtual void insertL(ut::KeyType key, int flavor, double* ptr) = 0;
            virtual void insertR(ut::KeyType key, int flavor, double* ptr) = 0;
            
            virtual ut::complex sign() const = 0;
            virtual void clean(itf::Hyb const&) = 0;
            
            virtual ~Bath() = default;
        };
        
    };
    

    template<typename HybVal>
    struct Bath : itf::Bath {
		Bath() : swapSign_(1), det_(1.), B_(0) {};
        Bath(Bath const& bath) = delete;
        Bath(Bath&& bath) = delete;
        Bath& operator=(Bath const& bath) = delete;
        Bath& operator=(Bath&& bath) = delete;
        ~Bath() = default;

		std::vector<Operator> const& opsL() const { return opsL_;};
		std::vector<Operator> const& opsR() const { return opsR_;};
        
        void insertL(ut::KeyType key, int flavor, double* ptr) {
            posL_[key] = opsL_.size(); opsL_.push_back(Operator(key, flavor, ptr));
        };
        void insertR(ut::KeyType key, int flavor, double* ptr) {
            posR_[key] = opsR_.size(); opsR_.push_back(Operator(key, flavor, ptr));
        };
		
        Matrix<HybVal>& B() { return B_;};
        Matrix<HybVal> const& B() const { return B_;};
        
        ut::complex sign() const { return det_/std::abs(det_);};
	
        void clean(itf::Hyb const& hybItf) {
            if(opsL_.size() != opsR_.size())
                throw std::runtime_error("bath::bath::clean: invalid configuration.");
            
            int const N = opsL_.size();
            
            if(N != B_.dim()) B_ = Matrix<HybVal>(N);
			  
            det_ = swapSign_;
        
			if(N) {
				Matrix<HybVal> toInvert(N);
				
                auto const& hyb = get<HybVal>(hybItf);
				for(int j = 0; j < N; ++j) 					
					for(int i = 0; i < N; ++i) 
						toInvert.at(i, j) = hyb(opsL_[i].flavor(), opsR_[j].flavor(), opsL_[i].key() - opsR_[j].key());
				
                std::fill(B_.data(), B_.data() + N*N, .0);
                for(int i = 0; i < N; ++i) B_.at(i, i) = 1.;
				
				int* ipiv = new int[N]; int info;
				gesv(&N, &N, toInvert.data(), &N, ipiv, B_.data(), &N, &info);
				for(int i = 0; i < N; ++i) 
					det_ *= (ipiv[i] != i + 1 ? -toInvert.at(i, i) : toInvert.at(i, i));
				delete[] ipiv;
			}
		}

	private:
		int swapSign_;
		HybVal det_;
		Matrix<HybVal> B_;
        
        std::vector<Operator> opsL_;
        std::vector<Operator> opsR_;
        std::map<ut::KeyType, int> posL_;
        std::map<ut::KeyType, int> posR_;
        
        template<typename, typename> friend struct Update;
        
        static double to_real(double arg) { return arg;};
        static double to_real(ut::complex const& arg) { return arg.real();};
	};
    
    
    template<typename HybVal> Bath<HybVal>& get(itf::Bath& bathItf) {
        return static_cast<Bath<HybVal>&>(bathItf);
    };
    
    template<typename HybVal> Bath<HybVal> const& get(itf::Bath const& bathItf) {
        return static_cast<Bath<HybVal> const&>(bathItf);
    };
}

#endif
