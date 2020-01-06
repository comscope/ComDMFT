#ifndef OBSERVABLES_DENSITYMATRIX_H
#define OBSERVABLES_DENSITYMATRIX_H

#include <vector>

#include "Observable.h"
#include "../Utilities.h"
#include "../Data.h"
#include "../State.h"
#include "../../../include/measurements/Measurements.h"

namespace obs {
    
    struct DMStatic {
        static std::string name() { return "density matrix"; };
        static double get(imp::itf::Dynamic const& dyn) { return 1.; };
    };
    
    struct DMDynamic {
        static std::string name() { return "density matrix dyn"; };
        static double get(imp::itf::Dynamic const& dyn) { return dyn.dyn(); };
    };
    
    template<typename Alloc, typename HybVal, typename Select>
    struct DensityMatrix : Observable {
        DensityMatrix() = delete;
        DensityMatrix(jsx::value const& jParams, data::Data const& data) :
        phase_(Phase::Sample),
        samples_(0) {
             init(imp::get<Alloc>(data.eig()), HybVal());
        };
        DensityMatrix(DensityMatrix const&) = delete;
        DensityMatrix(DensityMatrix&&) = delete;
        DensityMatrix& operator=(DensityMatrix const&) = delete;
        DensityMatrix& operator=(DensityMatrix&&) = delete;
        ~DensityMatrix() = default;
        
        bool sample(ut::complex const sign, data::Data const& data, state::State& state, imp::itf::Batcher& batcher) {
            switch(phase_) {
                case Phase::Sample:
                    for(auto sec : state.densityMatrix())
                        acc(sign, Select::get(state.dyn()), imp::get<Alloc>(state.densityMatrix()), sec, batcher, HybVal());
                
                    phase_ = Phase::Finalize; return false;
                    
                case Phase::Finalize:
                    phase_ = Phase::Sample;
            }
            
            ++samples_; return true;
        }
        
        void store(data::Data const&, jsx::value&, std::int64_t) {};
        
        void finalize(data::Data const& data, jsx::value& measurements, std::int64_t) {
            if(!measurements[std::string(Select::name())].is<jsx::array_t>())
                measurements[Select::name()] = jsx::array_t(data.eig().sectorNumber());
            
            for(int sec = 0; sec < data.eig().sectorNumber(); ++sec)
                store(imp::get<Alloc>(data.eig()), measurements[Select::name()], sec, HybVal());
        };

    private:
        enum class Phase { Sample, Finalize };
        
        Phase phase_;
        std::int64_t samples_;
        std::vector<std::unique_ptr<imp::Matrix<Alloc>>> accDensityMatrixR_;
        std::vector<std::unique_ptr<imp::Matrix<Alloc>>> accDensityMatrixI_;
        
        
        void init(imp::EigenValues<Alloc> const& eig, double) {
            accDensityMatrixR_.resize(eig.sectorNumber() + 1);
            for(int sec = eig.sectorNumber(); sec; --sec)
                accDensityMatrixR_[sec].reset(new imp::Matrix<Alloc>(typename imp::Matrix<Alloc>::Zero(eig.at(sec).dim())));
        };
        
        void acc(ut::complex const sign, double fact, imp::DensityMatrix<Alloc> const& densityMatrix, int sec, imp::itf::Batcher& batcher, double) {
            imp::axpy(*accDensityMatrixR_[sec], sign.real()*fact/densityMatrix.Z(), densityMatrix.mat(sec), batcher);
        };
        
        void store(imp::EigenValues<Alloc> const& eig, jsx::value& measurements, int sec, double) {
            int const dim = eig.at(sec + 1).dim();
            std::vector<double> temp(dim*dim, .0);
            imp::axpy(temp.data(), 1., *accDensityMatrixR_[sec + 1]);
            
            int const dim0 = eig.at(sec + 1).dim0();
            std::vector<double> temp0(dim0*dim0, .0);
            for(int i = 0; i < dim; ++i)
                for(int j = 0; j < dim; ++j)
                    temp0[i*dim0 + j] = (temp[i*dim + j] + temp[j*dim + i])/2.;
            
            measurements[sec] << meas::fix(temp0, samples_);
        };
        
        
        void init(imp::EigenValues<Alloc> const& eig, ut::complex) {
            accDensityMatrixR_.resize(eig.sectorNumber() + 1);
            for(int sec = eig.sectorNumber(); sec; --sec)
                accDensityMatrixR_[sec].reset(new imp::Matrix<Alloc>(typename imp::Matrix<Alloc>::Zero(eig.at(sec).dim())));
            
            accDensityMatrixI_.resize(eig.sectorNumber() + 1);
            for(int sec = eig.sectorNumber(); sec; --sec)
                accDensityMatrixI_[sec].reset(new imp::Matrix<Alloc>(typename imp::Matrix<Alloc>::Zero(eig.at(sec).dim())));
        };
        
        void acc(ut::complex const sign, double fact, imp::DensityMatrix<Alloc> const& densityMatrix, int sec, imp::itf::Batcher& batcher, ut::complex) {
            imp::axpy(*accDensityMatrixR_[sec], sign.real()*fact/densityMatrix.Z(), densityMatrix.mat(sec), batcher);
            imp::axpy(*accDensityMatrixI_[sec], sign.imag()*fact/densityMatrix.Z(), densityMatrix.mat(sec), batcher);
        };
        
        void store(imp::EigenValues<Alloc> const& eig, jsx::value& measurements, int sec, ut::complex) {
            int const dim = eig.at(sec + 1).dim();
            std::vector<double> tempR(dim*dim, .0);
            imp::axpy(tempR.data(), 1., *accDensityMatrixR_[sec + 1]);
            std::vector<double> tempI(dim*dim, .0);
            imp::axpy(tempI.data(), 1., *accDensityMatrixI_[sec + 1]);
            
            int const dim0 = eig.at(sec + 1).dim0();
            std::vector<ut::complex> temp0(dim0*dim0, .0);
            for(int i = 0; i < dim; ++i)
                for(int j = 0; j < dim; ++j)
                    temp0[i*dim0 + j] = ut::complex((tempR[i*dim + j] + tempR[j*dim + i])/2., (tempI[i*dim + j] - tempI[j*dim + i])/2.);
            
            measurements[sec] << meas::fix(temp0, samples_);
        };
    };
}

#endif
