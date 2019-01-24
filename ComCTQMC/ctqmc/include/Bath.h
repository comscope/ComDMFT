#ifndef BATH
#define BATH

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "Utilities.h"
#include "Updates.h"

#include "../../include/LinAlg.h"

//Vielleicht insert und erase noch verändern wie in CTBI hürütütütüt pluschter pluschter

namespace ba {
	struct Operator {
		Operator() = delete;
        Operator(ut::KeyType key, int flavor) : key_(key), flavor_(flavor) {};
        Operator(Operator const&) = delete;
        Operator(Operator&&) = default;
        Operator& operator=(Operator const&) = default;
        Operator& operator=(Operator&&) = delete;
        
        ut::KeyType& key() { return key_;};
        ut::KeyType key() const { return key_;};
        
        int flavor() { return flavor_;};
        int flavor() const { return flavor_;};
        
    private:
        ut::KeyType key_; int flavor_;
	};
	
	struct Matrix {
        Matrix() = delete;
		Matrix(int dim) : dim_(dim), data_(dim_*dim_) {};
        Matrix(Matrix const&) = delete;
        Matrix(Matrix&&) = delete;
        Matrix& operator=(Matrix const&) = delete;
        Matrix& operator=(Matrix&& other) = default;
		double& at(int i, int j) { return data_[i + dim_*j];};
		double const& at(int i, int j) const { return data_[i + dim_*j];};
		double* data() { return data_.data();};
		double const* data() const { return data_.data();};
		double* data(int i, int j) { return data_.data() + i + j*dim_;};
		double const* data(int i, int j) const { return data_.data() + i + j*dim_;};
        ~Matrix() = default;
	private:
        int dim_;
        std::vector<double> data_;
	};

	struct Bath { 
		Bath() : swapSign_(1), det_(1.), B_(0) {};
		template<class Hyb>
        Bath(jsx::object const& jConfig, Hyb const& hyb) : swapSign_(1), B_(jConfig.at("keyL").size()) {
			auto const& keyL = jConfig.at("keyL").array(); auto const& flavorL = jConfig.at("flavorL").array();
			auto const& keyR = jConfig.at("keyR").array(); auto const& flavorR = jConfig.at("flavorR").array();
			
			std::size_t const N = keyL.size();
			if(flavorL.size() != N || keyR.size() != N || flavorR.size() != N)
				throw std::runtime_error("Bath: missmatch in size of configurations.");
			
			for(std::size_t i = 0; i < N; ++i) {
				opsL_.push_back(Operator(keyL[i].int64(), flavorL[i].int64()));
				opsR_.push_back(Operator(keyR[i].int64(), flavorR[i].int64()));
				posL_[keyL[i].int64()] = posR_[keyR[i].int64()] = i;
			}
			
            clean(hyb);
		}
        Bath(Bath const& bath) = delete;
        Bath(Bath&& bath) = delete;
        Bath& operator=(Bath const& bath) = delete;
        Bath& operator=(Bath&& bath) = default;

		std::vector<Operator> const& opsL() const { return opsL_;};
		std::vector<Operator> const& opsR() const { return opsR_;};
		
        Matrix& B() { return B_;};
        Matrix const& B() const { return B_;};
        
		double det() const { return det_;};
	
		template<class Hyb>
		void clean(Hyb const& hyb) {
			int const N = opsL_.size();
			det_ = swapSign_;
			
			if(N) {
				Matrix toInvert(N);
				
				for(int j = 0; j < N; ++j) 					
					for(int i = 0; i < N; ++i) 
						toInvert.at(i,j) = hyb(opsL_[i].flavor(), opsR_[j].flavor(), opsL_[i].key() - opsR_[j].key());
				
				int const inc0 = 0; int const inc1 = 1;
				int const diagInc = N + 1; int const size = N*N;
				double const zero = .0; double const one = 1.;
				
				dcopy_(&size, &zero, &inc0, B_.data(), &inc1);
				dcopy_(&N, &one, &inc0, B_.data(), &diagInc);
				
				int* ipiv = new int[N]; int info;
				dgesv_(&N, &N, toInvert.data(), &N, ipiv, B_.data(), &N, &info);
				for(int i = 0; i < N; ++i) 
					det_ *= (ipiv[i] != i + 1 ? -toInvert.at(i, i) : toInvert.at(i, i));
				delete[] ipiv;
			}
		}
		
        void save(jsx::object& jConfig) const {
            jsx::array keyL(opsL_.size()); jsx::array flavorL(opsL_.size());
            jsx::array keyR(opsL_.size()); jsx::array flavorR(opsL_.size());
			
			for(std::size_t i = 0; i < opsL_.size(); ++i) {
                keyL[i] = opsL_[i].key(); flavorL[i] = static_cast<std::int64_t>(opsL_[i].flavor());
                keyR[i] = opsR_[i].key(); flavorR[i] = static_cast<std::int64_t>(opsR_[i].flavor());
			}
			
            jConfig["keyL"] = std::move(keyL); jConfig["flavorL"] = std::move(flavorL);
            jConfig["keyR"] = std::move(keyR); jConfig["flavorR"] = std::move(flavorR);
		};
		
		~Bath() = default;
	private:
		int swapSign_;
		double det_;
		Matrix B_;
        
        std::vector<Operator> opsL_;
        std::vector<Operator> opsR_;
        std::map<ut::KeyType, int> posL_;
        std::map<ut::KeyType, int> posR_;
        
        template<typename U> friend struct Updates;
	};
    
    template<typename> struct Updates {};
    
    template<>
    struct Updates<up::InsertTwo> {
        template<class Hyb>
        double ratio(up::InsertTwo const& u, std::vector<Bath> const& baths, Hyb const& hyb) {
            Bath const& bath = baths[u.bath];
            
            int const N = bath.opsL_.size();
            val_ = hyb(u.flavorL, u.flavorR, u.keyL() - u.keyR());
            
            if(N) {
                Bv_.resize(N); vec_.resize(N);
                for(int n = 0; n < N; ++n) vec_[n] = hyb(bath.opsL_[n].flavor(), u.flavorR, bath.opsL_[n].key() - u.keyR());
                
                char const no = 'n';
                int const inc = 1;
                double const zero = .0;
                double const one = 1.;
                dgemv_(&no, &N, &N, &one, bath.B_.data(), &N, vec_.data(), &inc, &zero, Bv_.data(), &inc);
                
                for(int n = 0; n < N; ++n) vec_[n] = hyb(u.flavorL, bath.opsR_[n].flavor(), u.keyL() - bath.opsR_[n].key());
                val_ -= ddot_(&N, vec_.data(), &inc, Bv_.data(), &inc);	
            }
            
            return val_;
        };

        void accept(up::InsertTwo const& u, std::vector<Bath>& baths) {
            Bath& bath = baths[u.bath];
            
            int const N = bath.opsL_.size();
            int const newN = N + 1;
            
            bath.posL_[u.keyL()] = bath.opsL_.size(); bath.opsL_.push_back(Operator(u.keyL(), u.flavorL));
            bath.posR_[u.keyR()] = bath.opsR_.size(); bath.opsR_.push_back(Operator(u.keyR(), u.flavorR));
            
            Matrix temp(newN);
            
            double fact = 1./val_;
            temp.at(N, N) = fact;
            
            if(N) {
                std::vector<double> hBTilde(N);
                
                char const yes = 't';
                int const inc = 1;
                double const zero = .0;
                double const one = 1.;
                dgemv_(&yes, &N, &N, &fact, bath.B_.data(), &N, vec_.data(), &inc, &zero, hBTilde.data(), &inc);
                dger_(&N, &N, &one, Bv_.data(), &inc, hBTilde.data(), &inc, bath.B_.data(), &N);
                
                for(int n = 0; n < N; ++n)
                    dcopy_(&N, bath.B_.data(0, n), &inc, temp.data(0, n), &inc);
                
                fact = -1./val_;
                dscal_(&N, &fact, Bv_.data(), &inc);
                dcopy_(&N, Bv_.data(), &inc, temp.data(0, N), &inc);
                
                double const minus = -1.;
                dscal_(&N, &minus, hBTilde.data(), &inc); 
                dcopy_(&N, hBTilde.data(), &inc, temp.data(N, 0), &newN); 
            }
            
            bath.B_ = std::move(temp);
            
            bath.det_ *= val_;
        };
        
        void reject(up::InsertTwo const& u, std::vector<Bath>& baths) {};
    private:
        std::vector<double> Bv_;
        std::vector<double> vec_;
        double val_;
    };

    template<>
    struct Updates<up::EraseTwo> {
        template<class Hyb>
        double ratio(up::EraseTwo const& u, std::vector<Bath> const& baths, Hyb const& hyb) {
            Bath const& bath = baths[u.bath];
            
            itL_ = bath.posL_.find(u.keyL()); itR_ = bath.posR_.find(u.keyR());
            
            return val_ = bath.B_.at(itR_->second, itL_->second);
        };
        
        void accept(up::EraseTwo const& u, std::vector<Bath>& baths) {
            Bath& bath = baths[u.bath];
            
            int const N = bath.opsL_.size(); int const newN = N - 1;
            int const posL = itL_->second; int const posR = itR_->second;
            
            Matrix temp(newN);
            
            if(newN) {
                int const inc = 1;
        
                if(posL != newN) {
                    dswap_(&N, bath.B_.data(0, newN), &inc, bath.B_.data(0, posL), &inc); bath.swapSign_ *= -1;
                    bath.posL_[bath.opsL_.back().key()] = posL; bath.opsL_[posL] = bath.opsL_.back();
                }
                if(posR != newN) {
                    dswap_(&N, bath.B_.data(newN, 0), &N, bath.B_.data(posR, 0), &N); bath.swapSign_ *= -1;
                    bath.posR_[bath.opsR_.back().key()] = posR; bath.opsR_[posR] = bath.opsR_.back();
                }
                
                for(int n = 0; n < newN; ++n)
                    dcopy_(&newN, bath.B_.data(0, n), &inc, temp.data(0, n), &inc);
                
                double const fact = -1./val_;
                dger_(&newN, &newN, &fact, bath.B_.data(0, newN), &inc, bath.B_.data(newN, 0), &N, temp.data(), &newN);
            }
            
            bath.B_ = std::move(temp);
            
            bath.posL_.erase(itL_); bath.opsL_.pop_back();
            bath.posR_.erase(itR_); bath.opsR_.pop_back();
            
            bath.det_ *= val_;
        };
        
        void reject(up::EraseTwo const& u, std::vector<Bath>& baths) {};
    private:
        std::vector<double> Bv_;
        std::vector<double> vec_;
        double val_;
        
        std::map<ut::KeyType, int>::const_iterator itL_;
        std::map<ut::KeyType, int>::const_iterator itR_;
    };
}

#endif
