#ifndef IMPURITY_OPERATORS_H
#define IMPURITY_OPERATORS_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <new>
#include <memory>


#include "Algebra.h"
#include "BitSet.h"
#include "Diagonal.h"
#include "../Utilities.h"
#include "../../../include/linalg/Operators.h"
#include "../../../include/mpi/Utilities.h"


namespace imp {
    
    namespace itf {
        
        struct Operator {
            virtual ~Operator() = default;
        };
        
        struct Operators {
            virtual int flavors() const = 0;
            virtual ~Operators() = default;
        };
        
    };
    
    //-----------------------------------------------------------------OPERATOR--------------------------------------------------------------------------------
    template<typename Alloc>
    struct Operator : itf::Operator {
        Operator() = delete;
        Operator(itf::EigenValues const& eigItf) :
        eig_(get<Alloc>(eigItf)),
        isMap_(eig_.sectorNumber() + 1),
        isMat_(eig_.sectorNumber() + 1),
        map_(new SectorNorm[eig_.sectorNumber() + 1]),
        mat_(static_cast<Matrix<Alloc>*>(::operator new(sizeof(Matrix<Alloc>)*(eig_.sectorNumber() + 1)))) {
        };
        Operator(char option, itf::EigenValues const& eig) : Operator(eig) {
            if(option != '1' && option != '0') throw std::runtime_error("Tr: option in operator constructor not defined");
            
            if(option == '1')
                for(int s = eig_.sectorNumber(); s; --s) {
                    set_map(s) = { s, .0 };
                    mat(s, typename Matrix<Alloc>::Identity(eig_.at(s).dim()));
                }
        };
        Operator(jsx::value& jOperator, itf::EigenValues const& eigItf) : Operator(eigItf) {
            if(static_cast<int>(jOperator.size()) != eig_.sectorNumber())
                throw(std::runtime_error("Tr: wrong number of sectors."));
            
            std::vector<int> temp(eig_.sectorNumber() + 1, 0); int start_sector = 1;
            for(auto& jBloc : jOperator.array()) {
                if(!jBloc("target").is<jsx::null_t>()) {
                    int target_sector = jBloc("target").int64() + 1;
                    
                    if(target_sector < 1 || eig_.sectorNumber() < target_sector)
                        throw std::runtime_error("Tr: target sector out of range.");
                    if(temp[target_sector]++)
                        throw std::runtime_error("Tr: target sector not unique.");
                    
                    auto const& matrix = jsx::at<io::rmat>(jBloc("matrix"));
                    
                    if(matrix.I() != eig_.at(target_sector).dim0() || matrix.J() != eig_.at(start_sector).dim0())
                        throw std::runtime_error("Tr: invalid matrix dimensions");
                    
                    set_map(start_sector) = { target_sector, .0 };
                    mat(start_sector, eig_.at(target_sector).dim(), eig_.at(start_sector).dim(), matrix);
                    
                    jBloc("matrix") = jsx::empty_t();
                } else
                    set_map(start_sector) = { 0, .0 };
                ++start_sector;
            }
        };
        Operator(Operator const&) = delete;
        Operator(Operator&&) = delete;
        Operator& operator=(Operator const&) = delete;
        Operator& operator=(Operator&&) = delete;
        ~Operator() {
            if(isMat_.any()) for(int s = eig_.sectorNumber(); s; --s) if(isMat_[s]) mat_[s].~Matrix();
            ::operator delete(mat_); delete [] map_;
        };
        
        int isMap(int s) const { return isMap_[s];};
        SectorNorm& set_map(int s) { isMap_.set(s); return map_[s];};
        
        SectorNorm& map(int s) {
            if(!isMap(s)) throw std::runtime_error("imp::Operator::map: invalid sector");
            return map_[s];
        };
        SectorNorm const& map(int s) const {
            if(!isMap(s)) throw std::runtime_error("imp::Operator::map const: invalid sector");
            return map_[s];
        };
        
        int isMat(int s) const { return isMat_[s];};
        template<typename... Args> Matrix<Alloc>& mat(int s, Args&& ... args) {
            if(isMat(s)) throw std::runtime_error("imp::Operator::mat: matrix allocated");
            new(mat_ + s) Matrix<Alloc>(std::forward<Args>(args)...);
            isMat_.set(s); return mat_[s];
        }
        
        Matrix<Alloc>& mat(int s) {
            if(!isMat_[s]) throw std::runtime_error("imp::Operator::mat: invalid sector");
            return mat_[s];
        };
        Matrix<Alloc> const& mat(int s) const {
            if(!isMat_[s]) throw std::runtime_error("imp::Operator::mat const: invalid sector");
            return mat_[s];
        };
        
        void assign(std::vector<int> const& sectors, SectorNormPtrs& arg) {
            arg.end = arg.begin;
            for(auto sec : sectors)
                if(!isMap(sec)) {
                    set_map(sec) = { sec, .0 };
                    *arg.end++ = map_ + sec;
                }
        };
        int missing(SectorNormPtrs& missing, SectorNormPtrs const& requested) {
            missing.begin = missing.end = requested.end;
            for(SectorNormPtrs::iterator it = requested.begin; it != requested.end; ++it)
                if(!isMap((*it)->sector)) {
                    set_map((*it)->sector) = { (*it)->sector, .0 }; // values are set by following map
                    *missing.end++ = map_ + (*it)->sector;
                }
            return missing.begin != missing.end;
        };
        void map(SectorNormPtrs& arg) const {
            SectorNormPtrs::iterator const end = arg.end; arg.end = arg.begin;
            for(SectorNormPtrs::iterator it = arg.begin; it != end; ++it) {
                (*it)->norm += map((*it)->sector).norm;
                if(((*it)->sector = map((*it)->sector).sector))
                    *arg.end++ = *it;
            }
        };
        
    private:
        EigenValues<Alloc> const& eig_;
        
        BitSet isMap_, isMat_;
        SectorNorm* const map_;
        Matrix<Alloc>* const mat_;
    };
    
    template<typename Alloc> Operator<Alloc>& get(itf::Operator& operatorItf) {
        return static_cast<Operator<Alloc>&>(operatorItf);
    };
    
    template<typename Alloc> Operator<Alloc> const& get(itf::Operator const& operatorItf) {
        return static_cast<Operator<Alloc> const&>(operatorItf);
    };
    
    
    
    template<typename Alloc>
    struct Operators : itf::Operators {
        Operators() = delete;
        Operators(jsx::value const& jParams, jsx::value& jOperators, itf::EigenValues const& eigItf) :
        flavors_(2*jOperators.size()),
        ops_(static_cast<Operator<Alloc>*>(::operator new(flavors_*sizeof(Operator<Alloc>)))) {
            mpi::cout << "Reading operators ... " << std::flush;
            
            int i = 0;
            for(auto& jOp : jOperators.array()) {
                jsx::value jOpDagg = linalg::transpose(jOp);
                
                new(ops_ + 2*i    ) Operator<Alloc>(jOp, eigItf);
                new(ops_ + 2*i + 1) Operator<Alloc>(jOpDagg, eigItf);
                
                ++i;
            }
            
            mpi::cout << "Ok" << std::endl;
        }
        Operators(Operators const&) = delete;
        Operators(Operators&&) = delete;
        Operators& operator=(Operators const&) = delete;
        Operators& operator=(Operators&&) = delete;
        ~Operators() {
            for(int f = 0; f < flavors_; ++f) ops_[f].~Operator();
            ::operator delete(ops_);
        };
        
        int flavors() const { return flavors_;};
        
        Operator<Alloc> const& at(int s) const { return ops_[s];};
        
    private:
        int const flavors_;
        Operator<Alloc>* ops_;
    };
    
    template<typename Alloc> Operators<Alloc>& get(itf::Operators& operatorsItf) {
        return static_cast<Operators<Alloc>&>(operatorsItf);
    };
    
    template<typename Alloc> Operators<Alloc> const& get(itf::Operators const& operatorsItf) {
        return static_cast<Operators<Alloc> const&>(operatorsItf);
    };
    
}

#endif  
