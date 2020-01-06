#ifndef OBSERVABLES_OPERATORS_H
#define OBSERVABLES_OPERATORS_H

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include "../Utilities.h"
#include "../impurity/Algebra.h"
#include "../impurity/Diagonal.h"
#include "../impurity/Operators.h"
#include "../../../include/JsonX.h"

namespace obs {
    
    namespace itf {
        
        struct BullaOperators {
            virtual int flavors() const = 0;
            virtual ~BullaOperators() = default;
        };
        
        struct Occupation {
            virtual int flavors() const = 0;
            virtual ~Occupation() = default;
        };
        
        struct BullaOccupation {
            virtual int flavors() const = 0;
            virtual ~BullaOccupation() = default;
        };
        
    };
    
    
    template<typename Alloc>
    struct BullaOperators : itf::BullaOperators {
        BullaOperators() = delete;
        BullaOperators(jsx::value& jInteraction, jsx::value& jOperators, imp::itf::EigenValues const& eigItf) :
        flavors_(2*jOperators.size()),
        ops_(static_cast<imp::Operator<Alloc>*>(::operator new(flavors_*sizeof(imp::Operator<Alloc>)))) {
            mpi::cout << "Reading bulla operators ... " << std::flush;

            if(static_cast<int>(jInteraction.size()) != eigItf.sectorNumber())
                throw(std::runtime_error("bu: wrong number of sectors in interaction."));
            
            int i = 0;
            for(auto& jOp : jOperators.array()) {
                jsx::value jBulla;
                
                linalg::mult('n', 'n',  1., jOp, jInteraction, .0, jBulla);
                linalg::mult('n', 'n', -1., jInteraction, jOp, 1., jBulla);
                
                jsx::value jBullaDagg = linalg::transpose(jBulla);
                
                new(ops_ + 2*i    ) imp::Operator<Alloc>(jBulla, eigItf);
                new(ops_ + 2*i + 1) imp::Operator<Alloc>(jBullaDagg, eigItf);
                
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
        imp::Operator<Alloc> const& at(int f) const {
            return ops_[f];
        };
        
    private:
        int const flavors_;
        imp::Operator<Alloc>* ops_;
    };
    
    template<typename Alloc> BullaOperators<Alloc>& get(itf::BullaOperators* bullaOpsItf) {
        return *static_cast<BullaOperators<Alloc>*>(bullaOpsItf);
    };
    
    template<typename Alloc> BullaOperators<Alloc> const& get(itf::BullaOperators const* bullaOpsItf) {
        return *static_cast<BullaOperators<Alloc> const*>(bullaOpsItf);
    };
    
    
    template<typename Alloc>
    struct Occupation : itf::Occupation {
        Occupation(jsx::value& jOperators, imp::itf::EigenValues const& eigItf) :
        flavors_(jOperators.size()),
        ops_(static_cast<imp::Operator<Alloc>*>(::operator new(flavors_*sizeof(imp::Operator<Alloc>)))) {
            mpi::cout << "Reading occupation ... " << std::flush;

            int i = 0;
            for(auto& jOp : jOperators.array()) {
                jsx::value jOcc;
                linalg::mult('t', 'n', 1., jOp, jOp, .0, jOcc);

                new(ops_ + i) imp::Operator<Alloc>(jOcc, eigItf);
                
                ++i;
            }
            
            mpi::cout << "Ok" << std::endl;
        };
        int flavors() const { return flavors_;};
        imp::Operator<Alloc> const& at(int f) const { return ops_[f];};
        ~Occupation() {
            for(int f = 0; f < flavors_; ++f) ops_[f].~Operator();
            ::operator delete(ops_);
        };
    private:
        int const flavors_;
        imp::Operator<Alloc>* ops_;
    };
    
    template<typename Alloc> Occupation<Alloc>& get(itf::Occupation* occItf) {
        return *static_cast<Occupation<Alloc>*>(occItf);
    };
    
    template<typename Alloc> Occupation<Alloc> const& get(itf::Occupation const* occItf) {
        return *static_cast<Occupation<Alloc> const*>(occItf);
    };
    
    
    template<typename Alloc>
    struct BullaOccupation : itf::BullaOccupation {
        BullaOccupation(jsx::value const& jParams, std::vector<double> const& filling, jsx::value jEigenValues, jsx::value& jOperators, imp::itf::EigenValues const& eigItf) :
        flavors_(jOperators.size()),
        ops_(static_cast<imp::Operator<Alloc>*>(::operator new(flavors_*sizeof(imp::Operator<Alloc>)))) {
            mpi::cout << "Reading bulla occupation ... " << std::flush;
            
            auto const mu = jParams("mu").real64();
            int sector = 1;
            
            for(auto& jEnergies : jEigenValues.array())
                for(auto& energy : jsx::at<io::rvec>(jEnergies))
                    energy += -mu*filling.at(sector);

            jsx::value jMatrixEigenValues = linalg::diag_to_operator(jEigenValues);
            
            int i = 0;
            for(auto& jOp : jOperators.array()) {
                jsx::value jOcc;
                linalg::mult('t', 'n', 1., jOp, jOp, .0, jOcc);
                
                jsx::value jBullaOcc;
                linalg::mult('n', 'n',  1., jMatrixEigenValues, jOcc, .0, jBullaOcc);
                linalg::mult('n', 'n', -1., jOcc, jMatrixEigenValues, 1., jBullaOcc);

                new(ops_ + i) imp::Operator<Alloc>(jBullaOcc, eigItf);
                
                ++i;
            }
            
            mpi::cout << "Ok" << std::endl;
        };
        int flavors() const { return flavors_;};
        imp::Operator<Alloc> const& at(int f) const { return ops_[f];};
        ~BullaOccupation() {
            for(int f = 0; f < flavors_; ++f) ops_[f].~Operator();
            ::operator delete(ops_);
        };
    private:
        int const flavors_;
        imp::Operator<Alloc>* ops_;
    };
    
    template<typename Alloc> BullaOccupation<Alloc>& get(itf::BullaOccupation* bullaOccItf) {
        return *static_cast<BullaOccupation<Alloc>*>(bullaOccItf);
    };
    
    template<typename Alloc> BullaOccupation<Alloc> const& get(itf::BullaOccupation const* bullaOccItf) {
        return *static_cast<BullaOccupation<Alloc> const*>(bullaOccItf);
    };
}

#endif
