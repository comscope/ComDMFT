#ifndef IMPURITY_DYNAMIC_H
#define IMPURITY_DYNAMIC_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>


#include "../Utilities.h"
#include "../../../include/mpi/Utilities.h"
#include "../../../include/JsonX.h"

//schä dömande a la lüünö si tü vulä ancor de mua e ell ma di wa tö fär encüle sche le gräk !! hahahhahahhahaha

namespace imp {
    
    namespace itf {
        
        struct Dynamic {
            Dynamic() = default;
            Dynamic(Dynamic const&) = delete;
            Dynamic(Dynamic&&) = delete;
            Dynamic& operator=(Dynamic const&) = delete;
            Dynamic& operator=(Dynamic&&) = delete;
            virtual ~Dynamic() = default;
            
            virtual void insert(ut::KeyType key, int type) {};
            virtual void erase(ut::KeyType key, int type) {};
            
            virtual ut::Zahl ratio() { return 1.;};
            
            virtual void accept() {};
            virtual void reject() {};
            
            virtual double E() const { return .0;};
            virtual double dyn() const { return .0;};
            
            virtual void clean() {};
        };
        
    };
    
    
    struct Simple {
        Simple() = delete;
        Simple(jsx::value const& jParams, jsx::value const& jDynF0) :
        nIt_(std::max(static_cast<int>((jParams.is("dynamic factor") ? jParams("dynamic factor").real64() : 4.)*jDynF0.size()), 1)),
        K_(nIt_ + 2),
        L_(nIt_ + 2) {
            mpi::cout << "Reading in dynamic interaction ... " << std::flush;  //scheisse dä kack durt en ewigkeit
            
            double const Ur0 = jDynF0(0).real64();

            for(int i = 0; i < nIt_ + 1; ++i) {
                double const tau = i/static_cast<double>(nIt_)*ut::beta();
                K_[i] = -.5*Ur0*tau*(ut::beta() - tau)/ut::beta();
                L_[i] = Ur0*tau/ut::beta();
            }
            for(std::size_t m  = 1; m < jDynF0.size(); ++m) {
                double const v_m = 2.*M_PI*m/ut::beta();
                double const F_m = jDynF0(m).real64();
                for(int i = 0; i < nIt_ + 1; ++i) {
                    double const tau = i/static_cast<double>(nIt_)*ut::beta();
                    K_[i] -= 2.*std::cos(v_m*tau)/(v_m*v_m*ut::beta())*F_m;
                    L_[i] += 2.*std::sin(v_m*tau)/(v_m*ut::beta())*F_m;
                }
            }
            K_.back() = .0;
            L_.back() = .0;
            
            mpi::cout << "Ok" << std::endl;
        };
        Simple(Simple const&) = delete;
        Simple(Simple&&) = delete;
        Simple& operator=(Simple const&) = delete;
        Simple& operator=(Simple&&) = delete;
        ~Simple() = default;

        double K(ut::KeyType key) const {
            double it = std::abs(key/static_cast<double>(ut::KeyMax))*nIt_; int i0 = static_cast<int>(it);
            return (1. - (it - i0))*K_[i0] + (it - i0)*K_[i0 + 1];
        };
        
        double L(ut::KeyType key) const {
            double const sign = key > 0 ? 1. : -1.;
            double it = std::abs(key/static_cast<double>(ut::KeyMax))*nIt_; int i0 = static_cast<int>(it);
            return sign*((1. - (it - i0))*L_[i0] + (it - i0)*L_[i0 + 1]);
        };
        
    private:
        int const nIt_;
        std::vector<double> K_;
        std::vector<double> L_;
    };

    
    struct Dynamic : itf::Dynamic {
        Dynamic() = delete;
        Dynamic(Simple const& func) :
        func_(func),
        w_(.0),
        wBackup_(.0),
        dyn_(nullptr) {
        };
        Dynamic(Dynamic const&) = delete;
        Dynamic(Dynamic&&) = delete;
        Dynamic& operator=(Dynamic const&) = delete;
        Dynamic& operator=(Dynamic&&) = delete;
        ~Dynamic() = default;
        
        void insert(ut::KeyType key, int type) { modOps_.push_back(Entry(key, 2*type - 1,  1));};
        void erase(ut::KeyType key, int type) { modOps_.push_back(Entry(key, 2*type - 1, -1));};
        
        ut::Zahl ratio() {
            for(std::vector<Entry>::iterator it = modOps_.begin(); it != modOps_.end(); ++it) {
                double temp = .0;
                
                for(std::vector<Entry>::iterator jt = ops_.begin(); jt != ops_.end(); ++jt)
                    temp += jt->sign*(func_.K(it->key - jt->key) + func_.K(jt->key - it->key));
                for(std::vector<Entry>::iterator jt = modOps_.begin(); jt != modOps_.end(); ++jt)
                    temp += jt->sign*jt->state*func_.K(it->key - jt->key);
                
                w_ += it->sign*it->state*temp;
            }
            
            return ut::exp((w_ - wBackup_)/2.);
        };
        
        void accept() {
            wBackup_ = w_;
            dyn_.reset(nullptr);
            
            for(std::vector<Entry>::iterator it = modOps_.begin(); it != modOps_.end(); ++it)
                if(it->state == 1)
                    ops_.push_back(*it);
                else {
                    std::vector<Entry>::iterator temp = std::find(ops_.begin(), ops_.end(), *it);
                    if(temp == ops_.end())
                        throw std::runtime_error("DynTr: Time not found");
                    ops_.erase(temp);
                }
            
            modOps_.clear();
        };
        void reject() {
            w_ = wBackup_;
            modOps_.clear();
        };
        
        double E() const { return w_/(2.*ut::beta());};
        
        double dyn() const {
            if(dyn_.get() == nullptr) {
                dyn_.reset(new double(.0));
                for(auto const& op : ops_) *dyn_ += op.sign*func_.L(op.key);
            }
            return *dyn_;
        };

        void clean() {
            w_ = .0;
            for(std::vector<Entry>::iterator it = ops_.begin(); it != ops_.end(); ++it)
                for(std::vector<Entry>::iterator jt = ops_.begin(); jt != ops_.end(); ++jt)
                    w_ += it->sign*jt->sign*func_.K(it->key - jt->key);
            wBackup_ = w_;
            dyn_.reset(nullptr);  // safer this way I think
        };
    
    private:
        struct Entry {
            Entry() {};
            Entry(ut::KeyType key, int sign, int state) : key(key), sign(sign), state(state) {};
            ut::KeyType key; int sign; int state;
            bool operator==(Entry const& rhs) { return key == rhs.key;};
        };
        
        Simple const& func_;
        
        double w_, wBackup_;
        mutable std::unique_ptr<double> dyn_;
        std::vector<Entry> modOps_, ops_;
    };
}

#endif
