#ifndef TRACEUTILITIES
#define TRACEUTILITIES

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <bitset>
#include <new>
#include <memory>
#include "Zahl.h"
#include "TraceAlgebra.h"


namespace tr {

	struct SectorNorm { int sector; double norm;};	
	struct SectorNormPtrs { typedef SectorNorm** iterator; SectorNorm** begin; SectorNorm** end;};
	
	struct ErrorBound { int sector; double norm; za::Zahl bound;};
	struct ErrorBounds { typedef ErrorBound* iterator; ErrorBound* begin; ErrorBound* end;};
	
	int compare(ErrorBound const& a, ErrorBound const& b) { return a.norm > b.norm;}
	
    // Koennte von std::vector abgeleited werden ... benötigt aber definition von move assignement (nothrow !) an verschiedenen stellen ... fuck that
    struct EigenValues  {
        EigenValues() = delete;
        EigenValues(jsx::value const& jParams, jsx::value& jEigenValues) :
        sectorNumber_(jEigenValues.size()),
        energies_(static_cast<Energies*>(::operator new(sizeof(Energies)*(sectorNumber_ + 1)))) {
            std::cout << "Reading in eigenvalues ... " << std::flush;
            int sector = 1;
            for(auto& jEnergies : jEigenValues.array())
                new(energies_ + sector++) Energies(jParams, jsx::at<io::rvec>(jEnergies));
            std::cout << "Ok" << std::endl;
        }
        EigenValues(EigenValues const&) = delete;
        EigenValues(EigenValues&&) = delete;
        EigenValues& operator=(EigenValues const&) = delete;
        EigenValues& operator=(EigenValues&&) = delete;
        
        int sectorNumber() const { return sectorNumber_;};
        Energies const& at(int s) const {
            return energies_[s];
        };
        
        ~EigenValues() {
            for(int s = sectorNumber_; s; --s) energies_[s].~Energies();
            ::operator delete(energies_);
        };
    private:
        int const sectorNumber_;
        Energies* energies_;
    };

	//----------------------------------------------------------------PROPAGATOR--------------------------------------------------------------------------

	struct Propagator {
        Propagator() = delete;
        Propagator(double time, EigenValues const& eig) : eig_(eig), time_(time), isProp_(eig_.sectorNumber() + 1), prop_(static_cast<Vector*>(::operator new(sizeof(Vector)*(eig_.sectorNumber() + 1)))) {}; /////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Propagator(Propagator const&) = delete;
        Propagator(Propagator&&) = delete;
        Propagator& operator=(Propagator const&) = delete;
        Propagator& operator=(Propagator&&) = delete;
        double time() const { return time_;};
		Vector const& at(int s) { 
			if(isProp_[s]) return prop_[s];
            new(prop_ + s) Vector(time_, eig_.at(s)); isProp_.flip(s); // Soetti ok sii so, oder ???
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
		EigenValues const& eig_;

		double const time_;
        ut::BitSet isProp_;     //eleganz vo arsch vo chue aber z'schnellschte won ich bis jetzt gfunde han
		Vector* const prop_;
	};
	
	//-------------------------------------------------------------------MATRIXBASE------------------------------------------------------------------------------
	struct BlocMatrixBase {
        int is(int s) { return isMat_[s];};

        template<typename... Args>
        Matrix& mat(int s, Args const & ... args) {
            new(mat_ + s) Matrix(args ...);
            isMat_.flip(s); return mat_[s];
        }
		Matrix& mat(int s) { return mat_[s];}; 
		Matrix const& mat(int s) const { return mat_[s];}; 
		
	protected:
        EigenValues const& eig_;
        
		SectorNorm* const secNorm_;
		
        BlocMatrixBase() = delete;
        BlocMatrixBase(EigenValues const& eig) :
        eig_(eig),
        secNorm_(new SectorNorm[eig_.sectorNumber() + 1]),
        isMat_(eig_.sectorNumber() + 1),
        mat_(static_cast<Matrix*>(::operator new(sizeof(Matrix)*(eig_.sectorNumber() + 1)))) {
        };
        BlocMatrixBase(BlocMatrixBase const&) = delete;
        BlocMatrixBase(BlocMatrixBase&&) = delete;
        BlocMatrixBase& operator=(BlocMatrixBase const&) = delete;
        BlocMatrixBase& operator=(BlocMatrixBase&&) = delete;
        
		~BlocMatrixBase() { 
			if(isMat_.any()) for(int s = eig_.sectorNumber(); s; --s) if(isMat_[s]) mat_[s].~Matrix();
            ::operator delete(mat_); delete [] secNorm_;
		};
    private:
        ut::BitSet isMat_;
        Matrix* const mat_;
	};

	
	//-----------------------------------------------------------------OPERATOR--------------------------------------------------------------------------------
	struct Operator : BlocMatrixBase {
        Operator() = delete;
        Operator(EigenValues const& eig) : BlocMatrixBase(eig), isSecNorm_(eig_.sectorNumber() + 1) {};
        Operator(char option, EigenValues const& eig) : BlocMatrixBase(eig), isSecNorm_(eig_.sectorNumber() + 1) {
            if(option != '1') throw std::runtime_error("Tr: option in operator constructor not defined");
            
			for(int s = eig_.sectorNumber(); s; --s) {
				sec(s) = s; norm(s) = .0; 
                mat(s, Matrix::Identity(eig.at(s).dim()));
			}
			sec(0) = 0;	
		};
        Operator(jsx::value& jOperator, EigenValues const& eig) : BlocMatrixBase(eig), isSecNorm_(eig_.sectorNumber() + 1) {
            if(static_cast<int>(jOperator.size()) != eig_.sectorNumber())
                throw(std::runtime_error("Tr: wrong number of sectors."));

            std::vector<int> temp(eig_.sectorNumber() + 1, 0); int start_sector = 1;
            for(auto& jBloc : jOperator.array()) {
                if(jBloc("target").type() != jsx::type::null) {
                    int target_sector = jBloc("target").int64() + 1;
                    
                    if(target_sector < 1 || eig_.sectorNumber() < target_sector)
                        throw std::runtime_error("Tr: target sector out of range.");
                    if(temp[target_sector]++)
                        throw std::runtime_error("Tr: target sector not unique.");
                    
                    auto const& matrix = jsx::at<io::rmat>(jBloc("matrix"));
                    
                    if(matrix.I() != eig.at(target_sector).dim0() || matrix.J() != eig.at(start_sector).dim0())
                        throw std::runtime_error("Tr: invalid matrix dimensions");
                    
                    sec(start_sector) = target_sector; norm(start_sector) = .0;
                    mat(start_sector, eig.at(target_sector).dim(), eig.at(start_sector).dim(), matrix);
                    
                    jBloc("matrix").reset();
                } else
                    sec(start_sector) = 0;
                ++start_sector;
            }
        };
        Operator(Operator const&) = delete;
        Operator(Operator&&) = delete;
        Operator& operator=(Operator const&) = delete;
        Operator& operator=(Operator&&) = delete;
        
		int missing(SectorNormPtrs& missing, SectorNormPtrs const& requested) {
			missing.begin = missing.end = requested.end;
			for(SectorNormPtrs::iterator it = requested.begin; it != requested.end; ++it)
				if(!isSecNorm_[(*it)->sector]) {
					isSecNorm_.flip((*it)->sector); 
					sec((*it)->sector) = (*it)->sector; norm((*it)->sector) = .0;
					 *missing.end++ = secNorm_ + (*it)->sector; 
				}
			return missing.begin != missing.end;	
		};
		
		void map(SectorNormPtrs& sectorNorms) const {
			SectorNormPtrs::iterator end = sectorNorms.end; sectorNorms.end = sectorNorms.begin;
			for(SectorNormPtrs::iterator it = sectorNorms.begin; it != end; ++it) {
				(*it)->norm += norm((*it)->sector); 
				if(((*it)->sector = sec((*it)->sector))) 
					*sectorNorms.end++ = *it; 
			}
		};
		
		int& sec(int s) { return secNorm_[s].sector;};
		int const& sec(int s) const { return secNorm_[s].sector;};

		double& norm(int s) { return secNorm_[s].norm;};
		double const& norm(int s) const { return secNorm_[s].norm;};
        
        ~Operator() = default;
	private:
        ut::BitSet isSecNorm_;
    };

    // Koennte von std::vector abgeleited werden ... benötigt aber definition von move assignement (nothrow !) an verschiedenen stellen ... fuck that
    struct Operators {
        Operators() = delete;
        Operators(jsx::value const& jParams, jsx::value& jOperators, EigenValues const& eig) :
        flavors_(2*jOperators.size()),
        identity_('1', eig),
        ops_(static_cast<Operator*>(::operator new(flavors_*sizeof(Operator)))) {
            std::cout << "Reading in operators ... ";  // put this back into data .....
            
            int i = 0;
            for(auto& jOp : jOperators.array()) {
                jsx::value jOpDagg; jOpDagg = jsx::array(jOp.size(), jsx::object{{"target", jsx::null()}});
                
                int start_sector = 0;
                for(auto& jBloc : jOp.array()) {
                    if(jBloc("target").type() != jsx::type::null) {
                        int target_sector = jBloc("target").int64();
                        
                        jOpDagg[target_sector]["target"] = static_cast<std::int64_t>(start_sector);
                        jOpDagg[target_sector]["matrix"] = jBloc("matrix");
                        
                        jsx::at<io::rmat>(jOpDagg[target_sector]["matrix"]).transpose();
                    }
                    ++start_sector;
                }
                
                new(ops_ + 2*i    ) Operator(jOp, eig);
                new(ops_ + 2*i + 1) Operator(jOpDagg, eig);
                
                ++i;
            }

            std::cout << "Ok" << std::endl;
        }
        Operators(Operators const&) = delete;
        Operators(Operators&&) = delete;
        Operators& operator=(Operators const&) = delete;
        Operators& operator=(Operators&&) = delete;
        
        Operator const& at(int s) const { return ops_[s];};
        Operator const& identity() const { return identity_;};
        
        ~Operators() {
            for(int f = 0; f < flavors_; ++f) ops_[f].~Operator();
            ::operator delete(ops_);
        };
    private:
        int const flavors_;        
        Operator identity_;
        Operator* ops_;
    };
    
	
	//-----------------------------------------------------------------DENSITYMATRIX---------------------------------------------------------------------------------
	struct DensityMatrix : BlocMatrixBase {
		typedef int const* iterator;
		
        DensityMatrix() = delete;
		DensityMatrix(SectorNormPtrs& ptrs, EigenValues const& eig) : BlocMatrixBase(eig), sec_(new int[eig_.sectorNumber() + 1]), end_(sec_), trace_(.0), sectorTrace_(new za::Zahl[eig_.sectorNumber() + 1]) {
			ptrs.end = ptrs.begin;
			for(int s = 1; s <= eig_.sectorNumber(); ++s) {
				*ptrs.end = secNorm_ + s; 
				(*ptrs.end)->sector = s;
				(*ptrs.end)->norm = .0;
				++ptrs.end;
			}
		};
        DensityMatrix(DensityMatrix const&) = delete;
        DensityMatrix(DensityMatrix&&) = delete;
        DensityMatrix& operator=(DensityMatrix const&) = delete;
        DensityMatrix& operator=(DensityMatrix&&) = delete;
		
		int surviving(ErrorBounds& errorBounds) {
			errorBounds.end = errorBounds.begin;
			for(int s = 1; s <= eig_.sectorNumber(); ++s)
				if(secNorm_[s].sector) {  
					errorBounds.end->sector = s; 
					errorBounds.end->norm = secNorm_[s].norm; 
					++errorBounds.end;
				}	
            return errorBounds.begin != errorBounds.end;
		};
        

        za::Zahl& Z(int s) {  *end_++ = s; return sectorTrace_[s];};
        za::Zahl& Z() { return trace_;};
        
		int const* begin() const { return sec_;};
		int const* end() const { return end_;};
        
		za::Zahl Z() const { return trace_;};
		double weight(int s) const { return (sectorTrace_[s]/trace_).toDouble();};
		
		~DensityMatrix() {
            delete[] sectorTrace_; delete[] sec_;
        };
	private:
		int* const sec_;
		int* end_;
		za::Zahl trace_;
		za::Zahl* const sectorTrace_;
	};
}

#endif  
