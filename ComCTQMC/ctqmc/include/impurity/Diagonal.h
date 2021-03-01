#ifndef CTQMC_INCLUDE_IMPURITY_DIAGONAL_H
#define CTQMC_INCLUDE_IMPURITY_DIAGONAL_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <new>
#include <memory>


#include "Algebra.h"
#include "BitSet.h"
#include "Dynamic.h"
#include "../Utilities.h"
#include "../../../include/mpi/Utilities.h"
#include "../../../include/JsonX.h"


namespace imp {
    
    namespace itf {
        
        struct EigenValues {
            virtual int sectorNumber() const = 0;
            virtual ~EigenValues() = default;
        };
        
    };
    
    
    struct SectorNorm { int sector; double norm;};
    struct SectorNormPtrs { typedef SectorNorm** iterator; SectorNorm** begin; SectorNorm** end;};
    
    // Koennte von std::vector abgeleited werden ... benötigt aber definition von move assignement (nothrow !) an verschiedenen stellen ... fuck that
    template<typename Mode>
    struct EigenValues : itf::EigenValues {
        EigenValues() = delete;
        EigenValues(jsx::value const& jParams, jsx::value jEigenValues, std::vector<double> const& filling, imp::Simple const* dynFunc) :
        sectorNumber_(jEigenValues.size()),
        energies_(static_cast<Energies<Mode>*>(::operator new(sizeof(Energies<Mode>)*(sectorNumber_ + 1)))) {
            mpi::cout << "Reading eigenvalues ... " << std::flush;
            
            auto const mu = jParams("mu").real64();
            
            int sector = 1;
            for(auto& jEnergies : jEigenValues.array()) {
                for(auto& energy : jsx::at<io::rvec>(jEnergies))
                    energy += -mu*filling.at(sector) + (dynFunc != nullptr ? dynFunc->shift(sector) : .0);
                new(energies_ + sector++) Energies<Mode>(jParams, jsx::at<io::rvec>(jEnergies));
            }
            
            mpi::cout << "Ok" << std::endl;
        }
        EigenValues(EigenValues const&) = delete;
        EigenValues(EigenValues&&) = delete;
        EigenValues& operator=(EigenValues const&) = delete;
        EigenValues& operator=(EigenValues&&) = delete;        
        ~EigenValues() {
            for(int s = sectorNumber_; s; --s) energies_[s].~Energies();
            ::operator delete(energies_);
        };
        
        int sectorNumber() const { return sectorNumber_;};
        Energies<Mode> const& at(int s) const {
            return energies_[s];
        };

    private:
        int const sectorNumber_;
        Energies<Mode>* energies_;
    };

	//----------------------------------------------------------------PROPAGATOR--------------------------------------------------------------------------
    template<typename Mode>
	struct Propagator {
        Propagator() = delete;
        Propagator(double time, EigenValues<Mode> const& eig) : eig_(eig), time_(time), isProp_(eig_.sectorNumber() + 1), prop_(static_cast<Vector<Mode>*>(::operator new(sizeof(Vector<Mode>)*(eig_.sectorNumber() + 1)))) {}; /////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Propagator(Propagator const&) = delete;
        Propagator(Propagator&&) = delete;
        Propagator& operator=(Propagator const&) = delete;
        Propagator& operator=(Propagator&&) = delete;
        double time() const { return time_;};
		Vector<Mode> const& at(int s) {
			if(isProp_[s]) return prop_[s];
            new(prop_ + s) Vector<Mode>(time_, eig_.at(s)); isProp_.set(s); // Soetti ok sii so, oder ???
            return prop_[s];
		};
        Vector<Mode> const& at(int s) const {
            if(!isProp_[s]) throw std::runtime_error("imp::Propagator::at: null pointer");
            return prop_[s];
        };
		void add(SectorNormPtrs& norms) const {
			for(SectorNormPtrs::iterator it = norms.begin; it != norms.end; ++it) 
			    (*it)->norm += time_*eig_.at((*it)->sector).min();
		};
		~Propagator() { 
			if(isProp_.any()) for(int s = eig_.sectorNumber(); s; --s) if(isProp_[s]) prop_[s].~Vector();
            ::operator delete(prop_);
		};
        
	private:
		EigenValues<Mode> const& eig_;

		double const time_;
        BitSet isProp_;     //eleganz vo arsch vo chue aber z'schnellschte won ich bis jetzt gfunde han
		Vector<Mode>* const prop_;
	};
    
    template<typename Mode> EigenValues<Mode>& get(itf::EigenValues& eigenValuesItf) {
        return static_cast<EigenValues<Mode>&>(eigenValuesItf);
    };
    
    template<typename Mode> EigenValues<Mode> const& get(itf::EigenValues const& eigenValuesItf) {
        return static_cast<EigenValues<Mode> const&>(eigenValuesItf);
    };
	
}

#endif  
