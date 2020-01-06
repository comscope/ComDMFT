#ifndef IMPURITY_DENSITYMATRIX_H
#define IMPURITY_DENSITYMATRIX_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <new>
#include <memory>

#include "Algebra.h"
#include "Diagonal.h"
#include "Operators.h"
#include "Product.h"
#include "../Utilities.h"
#include "../../../include/mpi/Utilities.h"


namespace imp {
    
    namespace itf {
        
        struct DensityMatrix {
            virtual ut::Flag surviving(itf::EigenValues const&) = 0;
            virtual ut::Flag decide(ut::Zahl const&, itf::Product&, itf::Batcher&) = 0;
            virtual std::vector<int>::const_iterator begin() const = 0;
            virtual std::vector<int>::const_iterator end() const = 0;
            virtual ut::Zahl Z() const = 0;
            virtual double weight(int) const = 0;
            virtual ~DensityMatrix() = default;
        };
        
    };
    
    struct Bound {
        int sec; double ln; ut::Zahl value;
    };
    
    inline int compare(Bound const& a, Bound const& b) {
        return a.ln > b.ln;
    };
    
    //-----------------------------------------------------------------DENSITYMATRIX---------------------------------------------------------------------------------
    template<typename Alloc>
    struct DensityMatrix : itf::DensityMatrix {
        typedef std::vector<int>::const_iterator iterator;
        
        DensityMatrix() = default;
        DensityMatrix(itf::Product& productItf, itf::EigenValues const& eigItf) : DensityMatrix(get<Alloc>(productItf), get<Alloc>(eigItf)) {
        };
        DensityMatrix(Product<Alloc>& product, EigenValues<Alloc> const& eig) :
        level_(product.height()),
        Z_(.0), z_(eig.sectorNumber() + 1) {
            std::vector<int> source(eig.sectorNumber()); std::iota(source.begin(), source.end(), 1);
            std::vector<int> target = product.map(product.first(), level_, source);

            for(int i = 0; i < eig.sectorNumber(); ++i)
                if(target[i]) {
                    if(target[i] == source[i])
                        bounds_.push_back({source[i], product.first()->op(level_)->map(source[i]).norm, ut::Zahl()});
                    else
                        throw std::runtime_error("imp::DensityMatrix::constructor: density matrix is not block-diagonal");
                }
        };
        DensityMatrix(DensityMatrix const&) = delete;
        DensityMatrix(DensityMatrix&&) = default;
        DensityMatrix& operator=(DensityMatrix const&) = delete;
        DensityMatrix& operator=(DensityMatrix&&) = default;
        ~DensityMatrix() = default;
        
        ut::Flag surviving(itf::EigenValues const& eig) {
            if(bounds_.size() == 0) return ut::Flag::Reject;

            for(auto it = bounds_.begin(); it != bounds_.end(); ++it) it->ln += get<Alloc>(eig).at(it->sec).ln_dim();

            std::sort(bounds_.begin(), bounds_.end(), &compare);
            
            for(auto it = bounds_.begin(); it != bounds_.end(); ++it) it->value = ut::exp(it->ln);
            for(auto it = bounds_.end() - 1; it != bounds_.begin(); --it) (it - 1)->value += it->value;
            bounds_.push_back({0, .0, ut::Zahl(.0)});
            
            bound_ = bounds_.begin(); return ut::Flag::Pending;
        };

        ut::Flag decide(ut::Zahl const& thresh, itf::Product& product, itf::Batcher& batcher) {
            if(ut::abs(Z_) + bound_->value <= ut::abs(thresh))
                return ut::Flag::Reject;
            else if(bound_->value <= std::numeric_limits<double>::epsilon()*ut::abs(Z_)) {
                op_ = get<Alloc>(product).first()->get_op(level_);
                return ut::Flag::Accept;
            }

            sectors_.push_back(bound_->sec);
            get<Alloc>(product).multiply(get<Alloc>(product).first(), level_, bound_->sec, batcher);
            trace(&z_[bound_->sec], &Z_, get<Alloc>(product).first()->op(level_)->mat(bound_->sec), batcher);
            
            ++bound_; return ut::Flag::Pending;
        };
        
        std::vector<int>::const_iterator begin() const { return sectors_.begin();};
        std::vector<int>::const_iterator end() const { return sectors_.end();};
        
        Matrix<Alloc> const& mat(int s) const { return static_cast<Operator<Alloc> const*>(op_.get())->mat(s);};
        
        ut::Zahl Z() const { return Z_;};
        double weight(int s) const { return (z_[s]/Z_).to_double();};
        
    private:
        std::unique_ptr<Operator<Alloc>> op_;
        
        int level_;
        std::vector<int> sectors_;
        std::vector<Bound> bounds_;
        std::vector<Bound>::iterator bound_;
        
        ut::Zahl Z_; std::vector<ut::Zahl> z_;
    };
    
    template<typename Alloc> DensityMatrix<Alloc>& get(itf::DensityMatrix& densityMatrixItf) {
        return static_cast<DensityMatrix<Alloc>&>(densityMatrixItf);
    };
    
    template<typename Alloc> DensityMatrix<Alloc> const& get(itf::DensityMatrix const& densityMatrixItf) {
        return static_cast<DensityMatrix<Alloc> const&>(densityMatrixItf);
    };
    
}

#endif  
