#ifndef DATA_H
#define DATA_H

#include <vector>

#include "Utilities.h"
#include "impurity/Diagonal.h"
#include "impurity/Operators.h"
#include "impurity/Dynamic.h"
#include "bath/Hyb.h"
#include "observables/Operators.h"

#include "../../include/mpi/Utilities.h"
#include "../../include/atomic/GenerateAtomic.h"
#include "../../include/options/Options.h"

namespace data {
    
    struct Data {
        Data() = delete;
        template<typename Alloc, typename HybVal>
        Data(jsx::value jParams, ut::Options<Alloc, HybVal>) {
            ut::beta = ut::Beta(jParams("beta").real64()); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            opt::complete_hloc(jParams);
            jParams("hloc") = ga::construct_hloc(jParams("hloc"));
            mpi::write(jParams("hloc"), "hloc.json");
            
            jParams["operators"] = ga::construct_annihilation_operators(jParams("hloc"));
            jParams("hybridisation")("functions") = mpi::read(jParams("hybridisation")("functions").string());
            if(jParams.is("dyn")) jParams("dyn") = mpi::read(jParams("dyn").string());
            
            filling_ = jsx::at<io::rvec>(jParams("hloc")("filling")); filling_.insert(filling_.begin(), 0);
            eig_.reset(new imp::EigenValues<Alloc>(jParams, filling(), jParams("hloc")("eigen values")));
            ide_.reset(new imp::Operator<Alloc>('1', eig()));
            
            if(jParams.is("green bulla") ? jParams("green bulla").boolean() : true) bullaOps_.reset(new obs::BullaOperators<Alloc>(jParams("hloc")("interaction"), jParams("operators"), eig()));
            if(jParams.is("occupation susceptibility direct") ? jParams("occupation susceptibility direct").boolean() : false) occ_.reset(new obs::Occupation<Alloc>(jParams("operators"), eig()));
            if(jParams.is("occupation susceptibility bulla")  ? jParams("occupation susceptibility bulla").boolean()  : false) bullaOcc_.reset(new obs::BullaOccupation<Alloc>(jParams, filling(), jParams("hloc")("eigen values"), jParams("operators"), eig()));
            ops_.reset(new imp::Operators<Alloc>(jParams, jParams("operators"), eig()));
            
            hyb_.reset(new bath::Hyb<HybVal>(jParams, jParams("hybridisation")("matrix"), jParams("hybridisation")("functions")));
            
            if(jParams.is("dyn")) dyn_.reset(new imp::Simple(jParams, jParams("dyn")));
            
        }
        Data(Data const&) = delete;
        Data(Data&&) = delete;
        Data& operator=(Data const&) = delete;
        Data& operator=(Data&&) = delete;
        ~Data() = default;

        std::vector<double> const& filling() const { return filling_;};
        imp::itf::EigenValues const& eig() const { return *eig_;}
        imp::itf::Operators const& ops() const { return *ops_;};
        imp::itf::Operator const& ide() const { return *ide_;};
        imp::Simple const* dyn() const { return dyn_.get();};
        
        bath::itf::Hyb const& hyb() const { return *hyb_;};
        
        obs::itf::BullaOperators const* bullaOps() const { return bullaOps_.get();};
        obs::itf::Occupation const* occ() const { return occ_.get();};
        obs::itf::BullaOccupation const* bullaOcc() const { return bullaOcc_.get();};

    private:
        std::vector<double> filling_;
        std::unique_ptr<imp::itf::EigenValues> eig_;
        std::unique_ptr<imp::itf::Operators> ops_;
        std::unique_ptr<imp::itf::Operator> ide_;
        std::unique_ptr<imp::Simple> dyn_;
        std::unique_ptr<bath::itf::Hyb> hyb_;
        std::unique_ptr<obs::itf::BullaOperators> bullaOps_;
        std::unique_ptr<obs::itf::Occupation> occ_;
        std::unique_ptr<obs::itf::BullaOccupation> bullaOcc_;
    };
}

#endif
