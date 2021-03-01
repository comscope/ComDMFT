#ifndef CTQMC_INCLUDE_IMPURITY_OPERATORS_H
#define CTQMC_INCLUDE_IMPURITY_OPERATORS_H

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
#include "../../../include/linalg/LinAlg.h"
#include "../../../include/linalg/Operators.h"
#include "../../../include/mpi/Utilities.h"


namespace imp {
    
    namespace itf {
        
        template<typename Value>
        struct Operator {
            virtual SectorNorm const& map(int) const = 0;
            virtual ~Operator() = default;
        };
        
        template<typename Value>
        struct Operators {
            virtual int flavors() const = 0;
            virtual Operator<Value> const& at(int) const = 0;
            virtual ~Operators() = default;
        };
        
    };
    
    //-----------------------------------------------------------------OPERATOR--------------------------------------------------------------------------------
    template<typename Mode, typename Value>
    struct Operator : itf::Operator<Value> {
        Operator() = delete;
        Operator(itf::EigenValues const& eigItf) :
        eig_(get<Mode>(eigItf)),
        isMap_(eig_.sectorNumber() + 1),
        isMat_(eig_.sectorNumber() + 1),
        map_(new SectorNorm[eig_.sectorNumber() + 1]),
        mat_(static_cast<Matrix<Mode, Value>*>(::operator new(sizeof(Matrix<Mode, Value>)*(eig_.sectorNumber() + 1)))) {
        };
        Operator(char const option, itf::EigenValues const& eigItf) : Operator(eigItf) {
            if(option == '1') {
                for(int s = eig_.sectorNumber(); s; --s) {
                    set_map(s) = { s, .0 };
                    mat(s, typename Matrix<Mode, Value>::Identity(eig_.at(s).dim()));
                }
            } else
                throw std::runtime_error("Tr: option in operator constructor not defined");
        };
        Operator(jsx::value const& jOperator, itf::EigenValues const& eigItf) : Operator(eigItf) {
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
                    
                    auto const& matrix = jsx::at<io::Matrix<Value>>(jBloc("matrix"));
                    
                    if(matrix.I() != eig_.at(target_sector).dim0() || matrix.J() != eig_.at(start_sector).dim0())
                        throw std::runtime_error("Tr: invalid matrix dimensions");
                    
                    double const norm = linalg::spectral_norm(matrix);
                    
                    if(norm != .0) {
                        set_map(start_sector) = { target_sector, std::log(norm) };
                        mat(start_sector, eig_.at(target_sector).dim(), eig_.at(start_sector).dim(), matrix);
                    } else
                        set_map(start_sector) = { 0, .0 };
                    
                    //jBloc("matrix") = jsx::empty_t();
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
        template<typename... Args> Matrix<Mode, Value>& mat(int s, Args&& ... args) {
            if(isMat(s)) throw std::runtime_error("imp::Operator::mat: matrix allocated");
            new(mat_ + s) Matrix<Mode, Value>(std::forward<Args>(args)...);
            isMat_.set(s); return mat_[s];
        }
        
        Matrix<Mode, Value>& mat(int s) {
            if(!isMat_[s]) throw std::runtime_error("imp::Operator::mat: invalid sector");
            return mat_[s];
        };
        Matrix<Mode, Value> const& mat(int s) const {
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
        EigenValues<Mode> const& eig_;
        
        BitSet isMap_, isMat_;
        SectorNorm* const map_;
        Matrix<Mode, Value>* const mat_;
    };
    
    template<typename Mode, typename Value> Operator<Mode, Value>& get(itf::Operator<Value>& operatorItf) {
        return static_cast<Operator<Mode, Value>&>(operatorItf);
    };
    
    template<typename Mode, typename Value> Operator<Mode, Value> const& get(itf::Operator<Value> const& operatorItf) {
        return static_cast<Operator<Mode, Value> const&>(operatorItf);
    };
    
    
    
    template<typename Mode, typename Value>
    struct Operators : itf::Operators<Value> {
        Operators() = delete;
        Operators(jsx::value const& jParams, jsx::value const& jOperators, itf::EigenValues const& eigItf) :
        flavors_(2*jOperators.size()),
        ops_(static_cast<Operator<Mode, Value>*>(::operator new(flavors_*sizeof(Operator<Mode, Value>)))) {
            mpi::cout << "Reading operators ... " << std::flush;
            
            int i = 0;
            for(auto& jOp : jOperators.array()) {
                jsx::value jOpDagg = linalg::conj<Value>(jOp);
                
                new(ops_ + 2*i    ) Operator<Mode, Value>(jOp, eigItf);
                new(ops_ + 2*i + 1) Operator<Mode, Value>(jOpDagg, eigItf);
                
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
        
        itf::Operator<Value> const& at(int i) const { return ops_[i];};
        
    private:
        int const flavors_;
        Operator<Mode, Value>* ops_;
    };
    
    template<typename Mode, typename Value> Operators<Mode, Value>& get(itf::Operators<Value>& operatorsItf) {
        return static_cast<Operators<Mode, Value>&>(operatorsItf);
    };
    
    template<typename Mode, typename Value> Operators<Mode, Value> const& get(itf::Operators<Value> const& operatorsItf) {
        return static_cast<Operators<Mode, Value> const&>(operatorsItf);
    };
    
}

#endif  
