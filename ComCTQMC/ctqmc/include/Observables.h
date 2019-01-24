#ifndef OBSERVABLES
#define OBSERVABLES

#include <vector>

#include "Utilities.h"
#include "Updates.h"
#include "Susc.h"
#include "Green.h"
#include "Data.h"
#include "Weight.h"

#include "../../include/Measurements.h"

namespace ob {
    
    struct DensityMatrix {
        DensityMatrix() = delete;
        DensityMatrix(tr::EigenValues const& eig) : accDensityMatrix_(eig.sectorNumber() + 1) {
            for(int s = eig.sectorNumber(); s; --s)
                accDensityMatrix_[s] = std::unique_ptr<tr::Matrix>(new tr::Matrix(tr::Matrix::Zero(eig.at(s).dim())));
        };
        DensityMatrix(DensityMatrix const&) = delete;
        DensityMatrix(DensityMatrix&&) = delete;
        DensityMatrix& operator=(DensityMatrix const&) = delete;
        DensityMatrix& operator=(DensityMatrix&&) = delete;
        
        void sample(double fact, tr::DensityMatrix const& densityMatrix) {
            for(auto sec : densityMatrix)
                tr::axpy(*accDensityMatrix_[sec], fact/densityMatrix.Z(), densityMatrix.mat(sec)); //Isch vo rechts nach links oder was. voll krassi sach.
        }
        
        void store(jsx::value& measurements, std::int64_t samples, tr::EigenValues const& eig) { //!!!!!!!!!!!!!!!!!!!!!! eig muess wo anders häre ...
            if(measurements.type() != jsx::type::array) measurements = jsx::array(eig.sectorNumber());
            
            for(int s = 0; s < eig.sectorNumber(); ++s) {
                int const dim = eig.at(s + 1).dim();
                std::vector<double> temp(dim*dim, .0);
                tr::axpy(temp.data(), 1., *accDensityMatrix_[s + 1]);
                
                int const dim0 = eig.at(s + 1).dim0();
                std::vector<double> temp0(dim0*dim0, .0);
                for(int i = 0; i < dim; ++i)
                    for(int j = 0; j < dim; ++j)
                        temp0[i*dim0 + j] = (temp[i*dim + j] + temp[j*dim + i])/2.;
                
                measurements[s] << meas::fix(temp0, samples);
            }
        }
        
        ~DensityMatrix() = default;
    private:
        std::vector<std::unique_ptr<tr::Matrix>> accDensityMatrix_;
    };
    
    
    template<class GreenMeas, class HybFunc>
    struct Observables {
        
        static Observables* Instance(jsx::value const& jParams, da::Data<HybFunc> const& data, we::Weight<HybFunc> const& weight) {
            if(counter_++ == 0) {
                if(jParams.is("density matrix") ? jParams("density matrix").boolean() : true) {
                    densityMatrix_ = new DensityMatrix(data.eig());
                    if(data.dyn() != nullptr) densityMatrixDyn_ = new DensityMatrix(data.eig());
                }
            }
            return new Observables(jParams, data, weight);
        };
        
        static void Destroy(Observables* instance) {
            delete instance;
            
            if(--counter_ == 0) {
                delete densityMatrixDyn_; densityMatrixDyn_ = nullptr;
                delete densityMatrix_; densityMatrix_ = nullptr;
            }
        };
        
    private:
        
        Observables() = delete;
        Observables(jsx::value const& jParams, da::Data<HybFunc> const& data, we::Weight<HybFunc> const& weight) :
        samples_(0),
        accSign_(.0),
        accE_(.0),
        acck_(.0),
        acckHist_((jParams.is("expansion histogram") ? jParams("expansion histogram").boolean() : true) ? 1 : 0),
        accDynDyn_(.0),
        green_(jParams),
        susc_(jParams, data.qns()) {
            for(auto const& x : weight.baths()) {
                auto const& opsL = x.opsL();
                auto const& opsR = x.opsR();
                
                for(std::size_t i = 0; i < opsL.size(); ++i) {
                    susc_.insert(opsR[i].key(), opsR[i].flavor());
                    susc_.insert(opsL[i].key(), opsL[i].flavor());
                }
            }
        };
        Observables(Observables const&) = delete;
        Observables(Observables&&) = delete;
        Observables& operator=(Observables const&) = delete;
        Observables& operator=(Observables&&) = delete;
        
        ~Observables() = default;
        
    public:
        
        void sample(da::Data<HybFunc> const& data, we::Weight<HybFunc> const& weight) {
            int const sign = weight.sign();
            
            accSign_ += sign;
            
            if(densityMatrixDyn_ != nullptr) {
                double temp = weight.dyn().dyn();
                accDynDyn_ += sign*temp*temp;
            }
            
            accE_ += sign*weight.dyn().E();
            for(tr::DensityMatrix::iterator it = weight.trace().densityMatrix().begin(); it != weight.trace().densityMatrix().end(); ++it) //!!!!!!!!!!!!!!!!!!!!
                tr::accE(&accE_, static_cast<double>(sign)/weight.trace().densityMatrix().Z(), weight.trace().densityMatrix().mat(*it), data.eig().at(*it));
            
            acck_ += sign*weight.trace().size()/2.;
            
            if(acckHist_.size()) {
                if(!(weight.trace().size()/2 < static_cast<int>(acckHist_.size())))
                    acckHist_.resize(weight.trace().size()/2 + 1, .0);
                acckHist_[weight.trace().size()/2] += sign;
            }
            
            susc_.sample(sign, weight.trace().densityMatrix());
            for(auto const& bath : weight.baths())
                green_.sample(sign, bath);
            
            ++samples_;
        };
        void store(jsx::value& measurements) {
            if(samples_) {
//                std::cout << "(Stream " + std::to_string(tr::Comm::iStream()) + ") k: " << acck_/samples_ << std::endl;
                
                measurements["Sign"] << meas::fix(accSign_, samples_); accSign_ = .0;
                measurements["Scal"]["E"] << meas::fix(accE_, samples_); accE_ = .0;
                measurements["Scal"]["k"] << meas::fix(acck_, samples_); acck_ = .0;
                
                if(acckHist_.size()) {
                    measurements["ExpansionHist"] << meas::var(acckHist_, samples_);
                    for(auto& x : acckHist_) x = .0;
                }
               
                if(densityMatrixDyn_ != nullptr) {
                    measurements["DynDyn"] << meas::fix(accDynDyn_, samples_); 
                    accDynDyn_ = .0;
                }

                susc_.store(measurements, samples_);
                green_.store(measurements, samples_);
            }
            
            samples_ = 0;
        };
        
        int csample(da::Data<HybFunc> const& data, we::Weight<HybFunc> const& weight) {
            int const sign = weight.sign();
            
            if(densityMatrix_ != nullptr) densityMatrix_->sample(sign, weight.trace().densityMatrix());
            if(densityMatrixDyn_ != nullptr) densityMatrixDyn_->sample(weight.dyn().dyn()*sign, weight.trace().densityMatrix());
            
            ++samplesCritical_;
            
            return 1;
        }
        void cstore(jsx::value& measurements, da::Data<HybFunc> const& data) {
            if(samplesCritical_) {
                if(densityMatrix_ != nullptr) densityMatrix_->store(measurements["DensityMatrix"], samplesCritical_, data.eig());
                if(densityMatrixDyn_ != nullptr) densityMatrixDyn_->store(measurements["DensityMatrixDyn"], samplesCritical_, data.eig());
            }
            
            samplesCritical_ = 0;
        };
        
        template<typename U>
        void update(U const& u) {
            up::get<U>(suscUpds_).update(susc_, u);
        };
        
        void clean() {
            susc_.clean();
        };
    private:
        static int counter_;
        
        static std::int64_t samplesCritical_;
        static DensityMatrix* densityMatrix_;          //delete
        static DensityMatrix* densityMatrixDyn_;       //delete
        
        std::int64_t samples_;
        double accSign_;
        
        double accE_, acck_;
        std::vector<double> acckHist_;
        double accDynDyn_;
        
        gr::Green<GreenMeas> green_;
        su::Susc susc_;
        
        up::Tuple<su::Updater, up::InsertTwo, up::EraseTwo> suscUpds_;
    };
    
    
    template<class GreenMeas, class HybFunc> int Observables<GreenMeas, HybFunc>::counter_ = 0;
    template<class GreenMeas, class HybFunc> std::int64_t Observables<GreenMeas, HybFunc>::samplesCritical_ = 0;
    template<class GreenMeas, class HybFunc> DensityMatrix* Observables<GreenMeas, HybFunc>::densityMatrix_ = nullptr;
    template<class GreenMeas, class HybFunc> DensityMatrix* Observables<GreenMeas, HybFunc>::densityMatrixDyn_ = nullptr;
}

#endif
