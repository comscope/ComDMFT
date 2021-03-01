#ifndef CTQMC_IMPURITY_OBSERVABLES_H
#define CTQMC_IMPURITY_OBSERVABLES_H

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include "Algebra.h"
#include "Diagonal.h"
#include "Operators.h"
#include "../../../include/JsonX.h"

namespace imp {
    
    namespace itf {
        
        template<typename Value>
        struct BullaOperators {
            virtual int flavors() const = 0;
            virtual ~BullaOperators() = default;
        };
        
        template<typename Value>
        struct Occupation {
            virtual int flavors() const = 0;
            virtual ~Occupation() = default;
        };
        
        template<typename Value>
        struct BullaOccupation {
            virtual int flavors() const = 0;
            virtual ~BullaOccupation() = default;
        };
        
    };
    
    
    template<typename Mode, typename Value>
    struct BullaOperators : itf::BullaOperators<Value> {
        BullaOperators() = delete;
        BullaOperators(jsx::value const& jInteraction, jsx::value const& jOperators, itf::EigenValues const& eig) :
        flavors_(2*jOperators.size()),
        ops_(static_cast<Operator<Mode, Value>*>(::operator new(flavors_*sizeof(Operator<Mode, Value>)))) {
            mpi::cout << "Reading bulla operators ... " << std::flush;
            
            if(static_cast<int>(jInteraction.size()) != eig.sectorNumber())
                throw(std::runtime_error("imp: wrong number of sectors in interaction."));
            
            int i = 0;
            for(auto& jOp : jOperators.array()) {
                jsx::value jBulla;
                
                linalg::mult<Value>('n', 'n',  1., jOp, jInteraction, .0, jBulla);
                linalg::mult<Value>('n', 'n', -1., jInteraction, jOp, 1., jBulla);
                
                jsx::value jBullaDagg = linalg::conj<Value>(jBulla);
                
                new(ops_ + 2*i    ) Operator<Mode, Value>(jBulla, eig);
                new(ops_ + 2*i + 1) Operator<Mode, Value>(jBullaDagg, eig);
                
                ++i;
            }
            
            mpi::cout << "Ok" << std::endl;
        };
        BullaOperators(BullaOperators const&) = delete;
        BullaOperators(BullaOperators&&) = delete;
        BullaOperators& operator=(BullaOperators const&) = delete;
        BullaOperators& operator=(BullaOperators&&) = delete;
        ~BullaOperators() {
            for(int f = 0; f < flavors_; ++f) ops_[f].~Operator();
            ::operator delete(ops_);
        };
        
        int flavors() const {
            return flavors_;
        };
        Operator<Mode, Value> const& at(int f) const {
            return ops_[f];
        };
        
    private:
        int const flavors_;
        Operator<Mode, Value>* ops_;
    };
    
    template<typename Mode, typename Value> BullaOperators<Mode, Value>& get(itf::BullaOperators<Value>& bullaOpsItf) {
        return static_cast<BullaOperators<Mode, Value>&>(bullaOpsItf);
    };
    
    template<typename Mode, typename Value> BullaOperators<Mode, Value> const& get(itf::BullaOperators<Value> const& bullaOpsItf) {
        return static_cast<BullaOperators<Mode, Value> const&>(bullaOpsItf);
    };
    
    
    template<typename Mode, typename Value>
    struct Occupation : itf::Occupation<Value> {
        Occupation() = delete;
        Occupation(jsx::value const& jOperators, itf::EigenValues const& eig) :
        flavors_(jOperators.size()),
        ops_(static_cast<Operator<Mode, Value>*>(::operator new(flavors_*sizeof(Operator<Mode, Value>)))) {
            mpi::cout << "Reading occupation ... " << std::flush;
            
            int i = 0;
            for(auto& jOp : jOperators.array()) {
                jsx::value jOcc;
                linalg::mult<Value>('c', 'n', 1., jOp, jOp, .0, jOcc);
                
                new(ops_ + i) Operator<Mode, Value>(jOcc, eig);
                
                ++i;
            }
            
            mpi::cout << "Ok" << std::endl;
        };
        Occupation(Occupation const&) = delete;
        Occupation(Occupation&&) = delete;
        Occupation& operator=(Occupation const&) = delete;
        Occupation& operator=(Occupation&&) = delete;
        ~Occupation() {
            for(int f = 0; f < flavors_; ++f) ops_[f].~Operator();
            ::operator delete(ops_);
        };
        
        int flavors() const { return flavors_;};
        Operator<Mode, Value> const& at(int f) const { return ops_[f];};
        
    private:
        int const flavors_;
        Operator<Mode, Value>* ops_;
    };
    
    template<typename Mode, typename Value> Occupation<Mode, Value>& get(itf::Occupation<Value>& occItf) {
        return static_cast<Occupation<Mode, Value>&>(occItf);
    };
    
    template<typename Mode, typename Value> Occupation<Mode, Value> const& get(itf::Occupation<Value> const& occItf) {
        return static_cast<Occupation<Mode, Value> const&>(occItf);
    };
    
    
    template<typename Mode, typename Value>
    struct BullaOccupation : itf::BullaOccupation<Value> {
        BullaOccupation() = delete;
        BullaOccupation(jsx::value const& jParams, std::vector<double> const& filling, jsx::value jEigenValues, jsx::value const& jOperators, itf::EigenValues const& eig) :
        flavors_(jOperators.size()),
        ops_(static_cast<Operator<Mode, Value>*>(::operator new(flavors_*sizeof(Operator<Mode, Value>)))) {
            mpi::cout << "Reading bulla occupation ... " << std::flush;
            
            auto const mu = jParams("mu").real64();
            int sector = 1;
            
            for(auto& jEnergies : jEigenValues.array())
                for(auto& energy : jsx::at<io::rvec>(jEnergies))
                    energy += -mu*filling.at(sector);
            
            jsx::value jMatrixEigenValues = linalg::diag_to_operator<Value>(jEigenValues);
            
            int i = 0;
            for(auto const& jOp : jOperators.array()) {
                jsx::value jOcc;
                linalg::mult<Value>('c', 'n', 1., jOp, jOp, .0, jOcc);
                
                jsx::value jBullaOcc;
                linalg::mult<Value>('n', 'n',  1., jMatrixEigenValues, jOcc, .0, jBullaOcc);
                linalg::mult<Value>('n', 'n', -1., jOcc, jMatrixEigenValues, 1., jBullaOcc);
                
                new(ops_ + i) Operator<Mode, Value>(jBullaOcc, eig);
                
                ++i;
            }
            
            mpi::cout << "Ok" << std::endl;
        };
        BullaOccupation(BullaOccupation const&) = delete;
        BullaOccupation(BullaOccupation&&) = delete;
        BullaOccupation& operator=(BullaOccupation const&) = delete;
        BullaOccupation& operator=(BullaOccupation&&) = delete;
        ~BullaOccupation() {
            for(int f = 0; f < flavors_; ++f) ops_[f].~Operator();
            ::operator delete(ops_);
        };
        
        int flavors() const { return flavors_;};
        Operator<Mode, Value> const& at(int f) const { return ops_[f];};
        
    private:
        int const flavors_;
        Operator<Mode, Value>* ops_;
    };
    
    template<typename Mode, typename Value> BullaOccupation<Mode, Value>& get(itf::BullaOccupation<Value>& bullaOccItf) {
        return static_cast<BullaOccupation<Mode, Value>&>(bullaOccItf);
    };
    
    template<typename Mode, typename Value> BullaOccupation<Mode, Value> const& get(itf::BullaOccupation<Value> const& bullaOccItf) {
        return static_cast<BullaOccupation<Mode, Value> const&>(bullaOccItf);
    };
    
}

#endif
