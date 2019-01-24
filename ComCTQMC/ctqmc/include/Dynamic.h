#ifndef DYNAMIC
#define DYNAMIC

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include "Utilities.h"
#include "Updates.h"

//schä dömande a la lüünö si tü vulä ancor de mua e ell ma di wa tö fär encüle sche le gräk !! hahahhahahhahaha

namespace dy {
    struct Simple {
        Simple() = delete;
        Simple(jsx::value const& jParams, jsx::value const& jDynF0) :
        nIt_(std::max(static_cast<int>((jParams.is("dynamic factor") ? jParams("dynamic factor").real64() : 4.)*jDynF0.size()), 1)),
        K_(nIt_ + 2),
        L_(nIt_ + 2) {
            std::cout << "Reading in dynamic interaction ... " << std::flush;  //scheisse dä kack durt en ewigkeit
            
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
            
            std::cout << "Ok" << std::endl;
        };
        Simple(Simple const&) = delete;
        Simple(Simple&&) = delete;
        Simple& operator=(Simple const&) = delete;
        Simple& operator=(Simple&&) = delete;
        
        double K(ut::KeyType key) const {
            double it = std::abs(key/static_cast<double>(ut::KeyMax))*nIt_; int i0 = static_cast<int>(it);
            return (1. - (it - i0))*K_[i0] + (it - i0)*K_[i0 + 1];
        };
        
        double L(ut::KeyType key) const {
            double it = key/static_cast<double>(ut::KeyMax)*nIt_; int i0 = static_cast<int>(it);
            return (1. - (it - i0))*L_[i0] + (it - i0)*L_[i0 + 1];
        };
        
        ~Simple() = default;
    private:
        int const nIt_;
        std::vector<double> K_;
        std::vector<double> L_;
    };
    
    
    struct Abstract {
        Abstract() {};
        
        virtual void insert(ut::KeyType key, int type) {};
        virtual void erase(ut::KeyType key, int type) {};
        
        virtual za::Zahl ratio() { return 1.;};
        
        virtual void accept() {};
        virtual void reject() {};
        
        virtual double E() const { return .0;};
        virtual double dyn() const { return .0;};
        
        virtual void clean() {};
        
        virtual ~Abstract() {};
    };
    
    struct Trace : Abstract {
        Trace() = delete;
        Trace(Simple const& func) :
        func_(func),
        w_(.0),
        wBackup_(.0) {
        };
        Trace(Trace const&) = delete;
        Trace(Trace&&) = delete;
        Trace& operator=(Trace const&) = delete;
        Trace& operator=(Trace&&) = delete;
        
        void insert(ut::KeyType key, int type) { modOps_.push_back(Entry(key, 2*type - 1,  1));};
        void erase(ut::KeyType key, int type) { modOps_.push_back(Entry(key, 2*type - 1, -1));};
        
        za::Zahl ratio() {
            for(std::vector<Entry>::iterator it = modOps_.begin(); it != modOps_.end(); ++it) {
                double temp = .0;
                
                for(std::vector<Entry>::iterator jt = ops_.begin(); jt != ops_.end(); ++jt)
                    temp += jt->sign*(func_.K(it->key - jt->key) + func_.K(jt->key - it->key));
                for(std::vector<Entry>::iterator jt = modOps_.begin(); jt != modOps_.end(); ++jt)
                    temp += jt->sign*jt->state*func_.K(it->key - jt->key);
                
                w_ += it->sign*it->state*temp;
            }
            
            return za::pow((w_ - wBackup_)/2.);
        };
        
        void accept() {
            wBackup_ = w_;
            
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
            double temp = .0;
            for(std::vector<Entry>::const_iterator it = ops_.begin(); it != ops_.end(); ++it) temp -= it->sign*func_.L(it->key);
            return temp;
        };

        void clean() {
            w_ = .0;
            for(std::vector<Entry>::iterator it = ops_.begin(); it != ops_.end(); ++it)
                for(std::vector<Entry>::iterator jt = ops_.begin(); jt != ops_.end(); ++jt)
                    w_ += it->sign*jt->sign*func_.K(it->key - jt->key);
            wBackup_ = w_;	
        };
        
        ~Trace() = default;
    private:
        struct Entry {
            Entry() {};
            Entry(ut::KeyType key, int sign, int state) : key(key), sign(sign), state(state) {};
            ut::KeyType key; int sign; int state;
            bool operator==(Entry const& rhs) { return key == rhs.key;};
        };
        
        Simple const& func_;
        
        double w_, wBackup_;
        
        std::vector<Entry> modOps_, ops_;
    };
    
    
    template<typename> struct Updates {};
    
    template<>
    struct Updates<up::InsertTwo> {
        za::Zahl ratio(up::InsertTwo const& u, Abstract& dyn) {
            dyn.insert(u.keyR(), u.flavorR%2);
            dyn.insert(u.keyL(), u.flavorL%2);
            return dyn.ratio();
        };
        
        void accept(Abstract& dyn) { dyn.accept();};
        void reject(Abstract& dyn) { dyn.reject();};
    };
    
    template<>
    struct Updates<up::EraseTwo> {
        za::Zahl ratio(up::EraseTwo const& u, Abstract& dyn) {
            dyn.erase(u.keyL(), u.flavorL%2);
            dyn.erase(u.keyR(), u.flavorR%2);
            return dyn.ratio();
        };
        
        void accept(Abstract& dyn) { dyn.accept();};
        void reject(Abstract& dyn) { dyn.reject();};
    };
}

#endif
