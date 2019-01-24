#ifndef SUSC
#define SUSC

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include "Utilities.h"
#include "Updates.h"

#include "../../include/Measurements.h"

namespace su {
    
    struct QuantumNumbers {
        QuantumNumbers() = delete;
        QuantumNumbers(jsx::value const& jParams, jsx::value& jQuantumNumbers, jsx::value const& jOperators, int const sectorNumber) :
        QNS_((sectorNumber + 1)*jQuantumNumbers.size()),
        qns_(2*jOperators.size()*jQuantumNumbers.size()) {
            for(auto& jQN : jQuantumNumbers.object()) {
                std::cout << "Checking quantum number " + jQN.first + " ... ";
                
                auto const& QN = jsx::at<io::rvec>(jQN.second);
                
                if(QN.size() != static_cast<std::size_t>(sectorNumber))
                    throw std::runtime_error("Susc: quantum numbers have wrong size.");
                
                for(std::size_t sector = 0; sector < QN.size(); ++sector)
                    QNS_[(sector + 1)*jQuantumNumbers.size() + names_.size()] = QN[sector];

                for(std::size_t i = 0; i < jOperators.size(); ++i) {
                    std::unique_ptr<double> qn; std::size_t start_sector = 0;
                    for(auto const& jBloc : jOperators(i).array()) {
                        if(jBloc("target").type() != jsx::type::null) {
                            std::size_t const target_sector = jBloc("target").int64();

                            if(qn.get() == nullptr)
                                qn = std::unique_ptr<double>(new double(QN[target_sector] - QN[start_sector]));
                            if(*qn != QN[target_sector] - QN[start_sector])
                                throw std::runtime_error("Susc: quantum numbers not supported !");
                        }
                        ++start_sector;
                    }
                    
                    qns_[(2*i    )*jQuantumNumbers.size() + names_.size()] =  *qn;
                    qns_[(2*i + 1)*jQuantumNumbers.size() + names_.size()] = -*qn;
                }
                
                names_.push_back(jQN.first);
                
                std::cout << "Ok" << std::endl;
            }
        };
        QuantumNumbers(QuantumNumbers const&) = delete;
        QuantumNumbers(QuantumNumbers&&) = delete;
        QuantumNumbers& operator=(QuantumNumbers const&) = delete;
        QuantumNumbers& operator=(QuantumNumbers&&) = delete;
        
        int size() const { return names_.size();};
        std::string const& names(int index) const { return names_[index];};
        double QN(int sector, int index) const { return QNS_[sector*names_.size() + index];};
        double qn(int flavor, int index) const { return qns_[flavor*names_.size() + index];};
        
        ~QuantumNumbers() = default;
    private:
        std::vector<std::string> names_;
        std::vector<double> QNS_;
        std::vector<double> qns_;
    };
    
    
    struct Susc {
        Susc() = delete;
        Susc(jsx::value const& jParams, QuantumNumbers const& quantumNumbers) :
        quantumNumbers_(quantumNumbers),
        nMat_(std::max(static_cast<int>(ut::beta()*jParams("susceptibility cutoff").real64()/(2*M_PI)), 1)),
        data_(.0, quantumNumbers_.size()*nMat_),
        accObs_(quantumNumbers_.size(), .0),
        accSusc_(quantumNumbers_.size(), std::vector<double>(nMat_, .0)) {
        };
        Susc(Susc const&) = delete;
        Susc(Susc&&) = delete;
        Susc& operator=(Susc const&) = delete;
        Susc& operator=(Susc&&) = delete;
        
        void insert(ut::KeyType key, int flavor) {
            Entry op(key, flavor);  update(op, 1);
            ops_.push_back(op);
        };
        void erase(ut::KeyType key, int flavor) {
            Entry op(key, flavor);  update(op, -1);
            std::vector<Entry>::iterator it = std::find(ops_.begin(), ops_.end(), op);
            if(it == ops_.end())
                throw std::runtime_error("Susc: Time not found");
            ops_.erase(it);
        };
        
        template<class D>
        void sample(int sign, D const& densityMatrix) {
            for(auto sector : densityMatrix) {
                double const fact = sign*densityMatrix.weight(sector)/ut::beta();
                for(int i = 0; i < quantumNumbers_.size(); ++i) {
                    double const data = ut::beta()*quantumNumbers_.QN(sector, i) - data_[i].real();
                    accObs_[i] += fact*data; accSusc_[i][0] += fact*data*data;
                }
            }
            
            for(std::size_t n = 1; n < nMat_; ++n) {
                double const fact = sign*ut::beta()/(4.*M_PI*M_PI*n*n);
                for(int i = 0; i < quantumNumbers_.size(); ++i)
                    accSusc_[i][n] += fact*(data_[n*quantumNumbers_.size() + i].real()*data_[n*quantumNumbers_.size() + i].real() + data_[n*quantumNumbers_.size() + i].imag()*data_[n*quantumNumbers_.size() + i].imag());
            }
        };
        void store(jsx::value& measurements, std::int64_t samples) {
            for(int i = 0; i < quantumNumbers_.size(); ++i) {
                measurements["Scal"][quantumNumbers_.names(i)] << meas::fix(accObs_[i], samples);  accObs_[i] = .0;
                measurements["Susc"][quantumNumbers_.names(i)] << meas::fix(accSusc_[i], samples);
                for(auto& x : accSusc_[i]) x = .0;
            }
        };
        
        void clean() {
            data_ = .0;
            for(std::vector<Entry>::const_iterator it = ops_.begin(); it != ops_.end(); ++it) update(*it, 1);
        };
        
        ~Susc() = default;
    private:
        struct Entry {
            Entry() {};
            Entry(ut::KeyType key, int flavor) : key(key), flavor(flavor) {};
            ut::KeyType key; int flavor;
            bool operator==(Entry const& rhs) { return key == rhs.key;};
        };
        
        QuantumNumbers const& quantumNumbers_;
        
        std::size_t const nMat_;
        std::valarray<ut::complex> data_;
        std::vector<double> accObs_;
        std::vector<std::vector<double> > accSusc_;
        std::vector<Entry> ops_;
        
        void update(Entry const& op, double state) {
            double const u = op.key/static_cast<double>(ut::KeyMax);
            
            for(int i = 0; i < quantumNumbers_.size(); ++i)
                data_[i] += state*ut::beta()*u*quantumNumbers_.qn(op.flavor, i);
            
            ut::complex const factExp = ut::complex(std::cos(2*M_PI*u), std::sin(2*M_PI*u)); ut::complex exp = 1.;
            for(std::size_t n = 1; n < nMat_; ++n) {
                exp *= factExp;
                
                ut::complex const fact = state*exp;
                for(int i = 0; i < quantumNumbers_.size(); ++i)
                    data_[quantumNumbers_.size()*n + i] += fact*quantumNumbers_.qn(op.flavor, i);
            }
        };
    };
    
    template<typename> struct Updater {};
    
    template<>
    struct Updater<up::InsertTwo> {
        void update(Susc& susc, up::InsertTwo const& u) {
            susc.insert(u.keyR(), u.flavorR);
            susc.insert(u.keyL(), u.flavorL);
        };
    };
    
    template<>
    struct Updater<up::EraseTwo> {
        void update(Susc& susc, up::EraseTwo const& u) {
            susc.erase(u.keyL(), u.flavorL);
            susc.erase(u.keyR(), u.flavorR);
        };
    };
}

#endif 
