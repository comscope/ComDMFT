#ifndef CTQMC_INCLUDE_IMPURITY_DENSITYMATRIX_H
#define CTQMC_INCLUDE_IMPURITY_DENSITYMATRIX_H

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
        
        template<typename Value>
        struct DensityMatrix {
            virtual ut::Flag surviving(itf::EigenValues const&) = 0;
            virtual ut::Flag decide(ut::Zahl<double> const&, itf::Product<Value>&, itf::Batcher<Value>&) = 0;
            virtual std::vector<int>::const_iterator begin() const = 0;
            virtual std::vector<int>::const_iterator end() const = 0;
            virtual ut::Zahl<Value> Z() const = 0;
            virtual Value weight(int) const = 0;
            virtual Value sign() const = 0;
            virtual ~DensityMatrix() = default;
        };
        
    };
    
    struct Bound {
        int sec; double ln; ut::Zahl<double> value;
    };
    
    inline int compare(Bound const& a, Bound const& b) {
        return a.ln > b.ln;
    };
    
    //-----------------------------------------------------------------DENSITYMATRIX---------------------------------------------------------------------------------
    template<typename Mode, typename Value>
    struct DensityMatrix : itf::DensityMatrix<Value> {
        typedef std::vector<int>::const_iterator iterator;
        
        DensityMatrix() = default;
        DensityMatrix(itf::Product<Value>& product, itf::EigenValues const& eig) :
        level_(product.height()),
        Z_(.0), z_(eig.sectorNumber() + 1) {
            std::vector<int> source(eig.sectorNumber()); std::iota(source.begin(), source.end(), 1);
            std::vector<int> target = get<Mode>(product).map(get<Mode>(product).first(), level_, source);

            for(int i = 0; i < eig.sectorNumber(); ++i)
                if(target[i] == source[i])
                    bounds_.push_back({source[i], get<Mode>(product).first()->op(level_)->map(source[i]).norm, ut::Zahl<double>()});
        };
        DensityMatrix(DensityMatrix const&) = delete;
        DensityMatrix(DensityMatrix&&) = default;
        DensityMatrix& operator=(DensityMatrix const&) = delete;
        DensityMatrix& operator=(DensityMatrix&&) = default;
        ~DensityMatrix() = default;
        
        ut::Flag surviving(itf::EigenValues const& eig) {
            if(bounds_.size() == 0) return ut::Flag::Reject;

            for(auto it = bounds_.begin(); it != bounds_.end(); ++it) it->ln += get<Mode>(eig).at(it->sec).ln_dim();

            std::sort(bounds_.begin(), bounds_.end(), &compare);
            
            for(auto it = bounds_.begin(); it != bounds_.end(); ++it) it->value = ut::exp(it->ln);
            for(auto it = bounds_.end() - 1; it != bounds_.begin(); --it) (it - 1)->value += it->value;
            bounds_.push_back({0, .0, ut::Zahl<double>(.0)});

            bound_ = bounds_.begin(); return ut::Flag::Pending;
        };

        ut::Flag decide(ut::Zahl<double> const& thresh, itf::Product<Value>& product, itf::Batcher<Value>& batcher) {
            if(ut::abs(Z_) + bound_->value <= ut::abs(thresh)) {
                return ut::Flag::Reject;
            } else if(bound_->value <= std::numeric_limits<double>::epsilon()*ut::abs(Z_)) {
                op_ = get<Mode>(product).first()->get_op(level_);
                return ut::Flag::Accept;
            }

            sectors_.push_back(bound_->sec);
            get<Mode>(product).multiply(get<Mode>(product).first(), level_, bound_->sec, batcher);
            trace(&z_[bound_->sec], &Z_, get<Mode>(product).first()->op(level_)->mat(bound_->sec), batcher);
            
            ++bound_; return ut::Flag::Pending;
        };
        
        std::vector<int>::const_iterator begin() const { return sectors_.begin();};
        std::vector<int>::const_iterator end() const { return sectors_.end();};
        
        Matrix<Mode, Value> const& mat(int s) const { return static_cast<Operator<Mode, Value> const*>(op_.get())->mat(s);};
        
        ut::Zahl<Value> Z() const { return Z_;};
        Value weight(int s) const { return (z_[s]/Z_).get();};
        
        Value sign() const { return Z_.mantissa()/std::abs(Z_.mantissa());};
        
    private:
        std::unique_ptr<Operator<Mode, Value>> op_;
        
        int level_;
        std::vector<int> sectors_;
        std::vector<Bound> bounds_;
        std::vector<Bound>::iterator bound_;
        
        ut::Zahl<Value> Z_; std::vector<ut::Zahl<Value>> z_;
    };
    
    template<typename Mode, typename Value> DensityMatrix<Mode, Value>& get(itf::DensityMatrix<Value>& densityMatrixItf) {
        return static_cast<DensityMatrix<Mode, Value>&>(densityMatrixItf);
    };
    
    template<typename Mode, typename Value> DensityMatrix<Mode, Value> const& get(itf::DensityMatrix<Value> const& densityMatrixItf) {
        return static_cast<DensityMatrix<Mode, Value> const&>(densityMatrixItf);
    };
    
}

#endif  
