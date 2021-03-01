#ifndef CTQMC_INCLUDE_IMPURITY_DYNAMIC_H
#define CTQMC_INCLUDE_IMPURITY_DYNAMIC_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>


#include "../Utilities.h"
#include "../../../include/atomic/Generate.h"
#include "../../../include/mpi/Utilities.h"
#include "../../../include/JsonX.h"

//schä dömande a la lüünö si tü vulä ancor de mua e ell ma di wa tö fär encüle sche le gräk !! hahahhahahhahaha


// Todo: get quantum number change as a double which should be given by the operator

namespace imp {
    
    namespace itf {
        
        struct Dynamic {
            Dynamic() = default;
            Dynamic(Dynamic const&) = delete;
            Dynamic(Dynamic&&) = delete;
            Dynamic& operator=(Dynamic const&) = delete;
            Dynamic& operator=(Dynamic&&) = delete;
            virtual ~Dynamic() = default;
            
            virtual void insert(ut::KeyType key, int flavor) {};
            virtual void erase(ut::KeyType key, int flavor) {};
            
            virtual ut::Zahl<double> ratio() { return 1.;};
            
            virtual void accept() {};
            virtual void reject() {};
            
            virtual double energy() const { return .0;};
            virtual double qkinks(int qn) const { return .0;};
            virtual double fkinks(int flavor, ut::KeyType key) const { return .0;};
            
            virtual void clean() {};
        };
        
    };
    
    
    //TODO: check if interaction matrix is symmetric (because the interaction is assumed to be real in this implementation) !
    
    
    struct Simple {
        Simple() = delete;
        Simple(jsx::value const& jParams, jsx::value jDyn) :
        size_(jDyn("quantum numbers").size()) {
            mpi::cout << "Reading in dynamic interaction ... " << std::flush;
            
            
            int const sectorNumber = jParams("hloc")("eigen values").size();
            int const flavors = jParams("operators").size();

            std::vector<std::vector<double>> q(size_);
            std::vector<std::vector<double>> Q(size_);
            std::vector<std::vector<double>> D0(size_, std::vector<double>(size_, .0));
            
            std::size_t nMat = std::numeric_limits<std::size_t>::max();
            
            if(jDyn("matrix").size() != size_)
                throw std::runtime_error("imp::Simple: matrix has wrong size !");
            
            for(int I = 0; I < size_; ++I)
            {
                q[I] = jsx::at<io::rvec>(jDyn("quantum numbers")(I));
                Q[I] = ga::construct_sector_qn(jParams("hloc"), jDyn("quantum numbers")(I));

                Q[I].insert(Q[I].begin(), 0);
                
                if(jDyn("matrix")(I).size() != size_)
                    throw std::runtime_error("imp::Simple: matrix has wrong size !");
                
                for(int J = 0; J < size_; ++J)
                    if(jDyn("matrix")(I)(J).string() != "")
                    {
                        auto const& D = jsx::at<io::rvec>(jDyn("functions")(jDyn("matrix")(I)(J).string()));
                        
                        nMat = std::min(nMat, D.size());  D0[I][J] = D[0];
                    }
            }
            
            if(nMat == std::numeric_limits<std::size_t>::max())
                throw std::runtime_error("imp::Simple: retarded interaction matrix is empty !");

            nIt_ = std::max(static_cast<int>((jParams.is("dynamic factor") ? jParams("dynamic factor").real64() : 4.)*nMat), 1);
            
            std::vector<std::vector<std::vector<double>>> Lqq(size_, std::vector<std::vector<double>>(size_, std::vector<double>(nIt_ + 2, .0)));
            std::vector<std::vector<std::vector<double>>> Kqq(size_, std::vector<std::vector<double>>(size_, std::vector<double>(nIt_ + 2, .0)));
            
            for(int I = 0; I < size_; ++I)
                for(int J = 0; J < size_; ++J)
                    if(jDyn("matrix")(I)(J).string() != "")
                    {
                        auto const& D = jsx::at<io::rvec>(jDyn("functions")(jDyn("matrix")(I)(J).string()));
                        
                        auto& L = Lqq[I][J];
                        auto& K = Kqq[I][J];
                        
                        for(int n = 0; n < nIt_ + 1; ++n)
                        {
                            double const tau = n/static_cast<double>(nIt_)*ut::beta();
                            
                            K[n] = -.5*D[0]*tau*(ut::beta() - tau)/ut::beta();
                            L[n] = D[0]*tau/ut::beta();
                            
                            for(std::size_t m  = 1; m < nMat; ++m)
                            {
                                double const v_m = 2.*M_PI*m/ut::beta();
                                
                                K[n] -= 2.*std::cos(v_m*tau)/(v_m*v_m*ut::beta())*D[m];
                                L[n] += 2.*std::sin(v_m*tau)/(v_m*ut::beta())*D[m];
                            }
                        }
                    }
            

            shift_.resize(sectorNumber + 1, .0);
            
            for(int I = 0; I < size_; ++I)
                for(int J = 0; J < size_; ++J)
                    for(int sec = 1; sec <= sectorNumber; ++sec)
                        shift_[sec] += .5*Q[I][sec]*D0[I][J]*Q[J][sec];
                        
            
            D0Qq_.resize(size_, std::vector<double>(sectorNumber + 1, .0));
            
            for(int I = 0; I < size_; ++I)
                for(int K = 0; K < size_; ++K)
                    for(int sec = 1; sec <= sectorNumber; ++sec)
                        D0Qq_[I][sec] += Q[K][sec]*D0[K][I];
            
            
            D0Qf_.resize(flavors, std::vector<double>(sectorNumber + 1, .0));
            
            for(int I = 0; I < size_; ++I)
                for(int i = 0; i < flavors; ++i)
                    for(int sec = 1; sec <= sectorNumber; ++sec)
                        D0Qf_[i][sec] += D0Qq_[I][sec]*q[I][i];
            

            Lfq_.resize(flavors, std::vector<std::vector<double>>(size_, std::vector<double>(nIt_ + 2, .0)));
            
            for(int I = 0; I < size_; ++I)
                for(int K = 0; K < size_; ++K)
                    for(int i = 0; i < flavors; ++i)
                        for(int n = 0; n < nIt_ + 1; ++n)
                            Lfq_[i][I][n] += q[K][i]*Lqq[K][I][n];

            
            Lff_.resize(flavors, std::vector<std::vector<double>>(flavors, std::vector<double>(nIt_ + 2, .0)));
            Kff_.resize(flavors, std::vector<std::vector<double>>(flavors, std::vector<double>(nIt_ + 2, .0)));
            
            for(int I = 0; I < size_; ++I)
                for(int J = 0; J < size_; ++J)
                    for(int i = 0; i < flavors; ++i)
                        for(int j = 0; j < flavors; ++j)
                            for(int n = 0; n < nIt_ + 1; ++n)
                            {
                                Lff_[i][j][n] += q[I][i]*Lqq[I][J][n]*q[J][j];
                                Kff_[i][j][n] += q[I][i]*Kqq[I][J][n]*q[J][j];
                            }

            mpi::cout << "Ok" << std::endl;
        };
        Simple(Simple const&) = delete;
        Simple(Simple&&) = delete;
        Simple& operator=(Simple const&) = delete;
        Simple& operator=(Simple&&) = delete;
        ~Simple() = default;
        
        int size() const {
            return size_;
        };
        
        double shift(int sec) const {
            return shift_[sec];
        };
        
        double D0Qq(int sec, int qn) const {
            return D0Qq_[qn][sec];
        };
        
        double D0Qf(int sec, int flavor) const {
            return D0Qf_[flavor/2][sec];
        };
        
        double Lfq(int flavor, int qn, ut::KeyType key) const {
            double it = std::abs(key/static_cast<double>(ut::KeyMax))*nIt_; int i0 = static_cast<int>(it);
            return (flavor%2 ? 1. : -1.)*(key > 0 ? 1. : -1.)*((1. - (it - i0))*Lfq_[flavor/2][qn][i0] + (it - i0)*Lfq_[flavor/2][qn][i0 + 1]);
        };
        
        double Lff(int flavorI, int flavorJ, ut::KeyType key) const {
            double it = std::abs(key/static_cast<double>(ut::KeyMax))*nIt_; int i0 = static_cast<int>(it);
            return (flavorI%2 ? 1. : -1.)*(key > 0 ? 1. : -1.)*((1. - (it - i0))*Lff_[flavorI/2][flavorJ/2][i0] + (it - i0)*Lff_[flavorI/2][flavorJ/2][i0 + 1]);
        };
        
        double Kff(int flavorI, int flavorJ, ut::KeyType key) const {
            double it = std::abs(key/static_cast<double>(ut::KeyMax))*nIt_; int i0 = static_cast<int>(it);
            return (flavorI%2 ? 1. : -1.)*(flavorJ%2 ? 1. : -1.)*((1. - (it - i0))*Kff_[flavorI/2][flavorJ/2][i0] + (it - i0)*Kff_[flavorI/2][flavorJ/2][i0 + 1]);
        };

    private:
        int const size_;
        
        int nIt_;
        std::vector<double> shift_;
        std::vector<std::vector<double>> D0Qq_;
        std::vector<std::vector<double>> D0Qf_;
        std::vector<std::vector<std::vector<double>>> Lfq_;
        std::vector<std::vector<std::vector<double>>> Lff_;
        std::vector<std::vector<std::vector<double>>> Kff_;
    };

    
    struct Dynamic : itf::Dynamic {
        Dynamic() = delete;
        Dynamic(Simple const& func) :
        func_(func),
        w_(.0),
        wBackup_(.0),
        qkinks_(func.size()) {
        };
        Dynamic(Dynamic const&) = delete;
        Dynamic(Dynamic&&) = delete;
        Dynamic& operator=(Dynamic const&) = delete;
        Dynamic& operator=(Dynamic&&) = delete;
        ~Dynamic() = default;
        
        void insert(ut::KeyType key, int flavor) {
            modOps_.push_back(Entry(key, flavor,  1));
        };
        void erase(ut::KeyType key, int flavor) {
            modOps_.push_back(Entry(key, flavor, -1));
        };
        
        ut::Zahl<double> ratio() {
            for(std::vector<Entry>::iterator it = modOps_.begin(); it != modOps_.end(); ++it) {
                double temp = .0;
                
                for(std::vector<Entry>::iterator jt = ops_.begin(); jt != ops_.end(); ++jt)
                    temp += func_.Kff(it->flavor, jt->flavor, it->key - jt->key) + func_.Kff(jt->flavor, it->flavor, jt->key - it->key);
                
                for(std::vector<Entry>::iterator jt = modOps_.begin(); jt != modOps_.end(); ++jt)
                    temp += jt->state*func_.Kff(it->flavor, jt->flavor, it->key - jt->key);
                
                w_ += it->state*temp;
            }
            
            return ut::exp((w_ - wBackup_)/2.);
        };
        
        void accept() {
            wBackup_ = w_;
            for(auto& qkink : qkinks_) qkink.reset(nullptr);
            
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
        
        double energy() const {
            return w_/(2.*ut::beta());
        };
        
        double qkinks(int qn) const {
            if(qkinks_[qn].get() == nullptr) {
                qkinks_[qn].reset(new double(.0));
                for(auto const& op : ops_) *qkinks_[qn] += func_.Lfq(op.flavor, qn, op.key);
            }
            return *qkinks_[qn];
        };
        
        double fkinks(int flavor, ut::KeyType key) const {
            double result = .0;
            for(auto const& op : ops_) result += func_.Lff(op.flavor, flavor, op.key - key);
            return result;
        };

        void clean() {
            w_ = .0;
            for(std::vector<Entry>::iterator it = ops_.begin(); it != ops_.end(); ++it)
                for(std::vector<Entry>::iterator jt = ops_.begin(); jt != ops_.end(); ++jt)
                    w_ += func_.Kff(it->flavor, jt->flavor, it->key - jt->key);
            wBackup_ = w_;
            for(auto& qkink : qkinks_) qkink.reset(nullptr);
        };
    
    private:
        struct Entry {
            Entry() {};
            Entry(ut::KeyType key, int flavor, int state) : key(key), flavor(flavor), state(state) {};
            ut::KeyType key; int flavor; int state;
            bool operator==(Entry const& rhs) { return key == rhs.key;};
        };
        
        Simple const& func_;
        
        double w_, wBackup_;
        std::vector<Entry> modOps_, ops_;
        
        mutable std::vector<std::unique_ptr<double>> qkinks_;
    };
}

#endif
