#ifndef DATA
#define DATA

#include <vector>

#include "Utilities.h"
#include "Trace.h"
#include "Dynamic.h"
#include "Susc.h"
#include "Hyb.h"

#include "../../include/MPIUtilities.h"
#include "../../include/impurity/GenerateAtomic.h"
#include "../../include/impurity/Hloc.h"
#include "../../include/impurity/Options.h"

namespace da {
    
    // No nit alles optimal im constructor .....
    
    template<typename HybFunc>
    struct Data {
        static Data const* Instance(jsx::value const& jParams) {
            if(counter_++ == 0) {
                //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ut::beta = ut::Beta(jParams("beta").real64());
                //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                jsx::value jImpurityModel;
                
                if(jParams.is("impurity"))
                    mpi::read(jParams("impurity").string(), jImpurityModel["impurity"]);
                else {
                    std::unique_ptr<Ga::GenerateAtomic> generateAtomic;

                    if(mpi::rank() == mpi::master) {
                        jImpurityModel["impurity"]["hloc"]["one body"] = jParams("hloc")("one body");
                        generateAtomic = std::unique_ptr<Ga::GenerateAtomic>(new Ga::GenerateAtomic(Ga::Hloc(options::get_hloc(jParams)), jImpurityModel("impurity")("hloc")));
                        mpi::broad_cast(jImpurityModel("impurity")("hloc"));
                    } else {
                        mpi::broad_cast(jImpurityModel["impurity"]["hloc"]);
                        generateAtomic = std::unique_ptr<Ga::GenerateAtomic>(new Ga::GenerateAtomic(jImpurityModel("impurity")("hloc")));
                    }
                    mpi::write(jImpurityModel("impurity")("hloc"), "Hloc.json");

                    generateAtomic->write_operators(jImpurityModel["impurity"]["operators"]);
                };
                
                jImpurityModel["hybridisation"] = jParams("hybridisation");
                mpi::read(jParams("hybridisation")("functions").string(), jImpurityModel["hybridisation"]["functions"]);
                
                if(jParams.is("dyn")) mpi::read(jParams("dyn").string(), jImpurityModel["dyn"]);
                
                auto const mu = jParams("mu").real64();
                auto const Ur0 = jImpurityModel.is("dyn") ? jImpurityModel("dyn")(0).real64() : .0;
                auto& jEigenValues = jImpurityModel("impurity")("hloc")("eigen values");
                auto& jN = jImpurityModel("impurity")("hloc")("quantum numbers")("N");
                
                if(jEigenValues.size() != jsx::at<io::rvec>(jN).size())
                    throw std::runtime_error("Data: missmatch in number of sectors in eigen values and quantum number N");
                
                for(std::size_t sector = 0; sector < jEigenValues.size(); ++sector)
                    for(auto& energy : jsx::at<io::rvec>(jEigenValues(sector)))
                        energy += -mu*jsx::at<io::rvec>(jN)[sector] + .5*Ur0*jsx::at<io::rvec>(jN)[sector]*jsx::at<io::rvec>(jN)[sector];
                
                data_ = new Data<HybFunc>(jParams, jImpurityModel);
            }
            
            ++counter_; return data_;
        }
        
        static void Destroy(Data const* data) {
            if(--counter_ == 0) {
                delete data_; data_ = nullptr;
            }
        }
        
    private:
        
        /// Idealerwiis soetti da jParams gar nit brucht werde ....
        Data() = delete;
        Data(jsx::value const& jParams, jsx::value& jImpurityModel) :
        eig_(jParams, jImpurityModel("impurity")("hloc")("eigen values")),
        ops_(jParams, jImpurityModel("impurity")("operators"), eig_),
        hyb_(jParams, jImpurityModel("hybridisation")("matrix"), jImpurityModel("hybridisation")("functions")),
        qns_(jParams, jImpurityModel("impurity")("hloc")("quantum numbers"), jImpurityModel("impurity")("operators"), eig_.sectorNumber()),  //scheisse das mit sectorNumber daa ... arsch vo chue.
        dyn_(jImpurityModel.is("dyn") ? new dy::Simple(jParams, jImpurityModel("dyn")) : nullptr) {
        }
        Data(Data const&) = delete;
        Data(Data&&) = delete;
        Data& operator=(Data const&) = delete;
        Data& operator=(Data&&) = delete;
        
        ~Data() = default;
        
    public:
        
        tr::EigenValues const& eig() const { return eig_;};
        tr::Operators const& ops() const { return ops_;};
        hy::Hyb<HybFunc> const& hyb() const { return hyb_;};
        su::QuantumNumbers const& qns() const { return qns_;};
        dy::Simple const* dyn() const { return dyn_.get();};
        
    private:
        static int counter_;
        static Data* data_;        //delete
        
        tr::EigenValues eig_;
        tr::Operators ops_;
        hy::Hyb<HybFunc> hyb_;
        su::QuantumNumbers qns_;
        std::unique_ptr<dy::Simple> dyn_;
    };
    
    template<class HybFunc> int Data<HybFunc>::counter_ = 0;
    template<class HybFunc> Data<HybFunc>* Data<HybFunc>::data_ = 0;
    
}

#endif
