#ifndef BATH_HYB_H
#define BATH_HYB_H

#include <cmath>
#include <iostream>
#include <sstream>

#include "../Utilities.h"
#include "../../../include/atomic/GenerateAtomic.h"

namespace bath {
    
    struct Block {
        std::vector<int>& flavorsL() { return flavorsL_;};
        std::vector<int>& flavorsR() { return flavorsR_;};
        
        std::vector<int> const& flavorsL() const { return flavorsL_;};
        std::vector<int> const& flavorsR() const { return flavorsR_;};
        
    private:
        std::vector<int> flavorsL_;
        std::vector<int> flavorsR_;
    };
    
    
    namespace itf {
        
        struct Hyb {
            virtual std::vector<Block> const& blocks() const = 0;
            virtual ~Hyb() = default;
        };
        
    };
    

    template<typename HybVal>
    struct Fit {
        Fit(double const beta, std::vector<std::complex<double>> const& hyb, std::vector<std::complex<double>> const& hybTransp) {
            if(hyb.size() != hybTransp.size()) throw std::runtime_error("bath::Fit: size if conjugate hybridisation functions not the same !");
            if(!hyb.size()) throw std::runtime_error("bath::Fit: no entries in hybridisation function !");
            
            N_ = std::max(static_cast<int>(hyb.size()*.1), 1);
            ut::complex D = .0, iwD = 0, D2 = .0, iwD2 = .0;
            for(std::size_t m = hyb.size() - N_; m < hyb.size(); ++m) {
                check(hyb[m], hybTransp[m], HybVal());
                
                ut::complex const iw{.0, M_PI*(2*m + 1)/beta};
                D    += hyb[m] + std::conj(hybTransp[m]);
                iwD  += iw*(hyb[m] - std::conj(hybTransp[m]));
                D2   += std::abs(hyb[m])*std::abs(hyb[m]) + std::abs(hybTransp[m])*std::abs(hybTransp[m]);
                iwD2 += iw*(std::abs(hyb[m])*std::abs(hyb[m]) - std::abs(hybTransp[m])*std::abs(hybTransp[m]));
            }
            
            moment_ = (iwD*D2 - D*iwD2)/(2.*N_*D2 - std::abs(D)*std::abs(D));
            eps_ = (2.*N_*iwD2 - iwD*std::conj(D))/(2.*N_*D2 - std::abs(D)*std::abs(D));
        }
        
        std::size_t N() const { return N_;};
        std::complex<double> moment() const { return moment_;};
        std::complex<double> eps() const { return eps_;};
    private:
        std::size_t N_;
        std::complex<double> moment_;
        std::complex<double> eps_;
        
        void check(ut::complex, ut::complex, ut::complex) {};
        void check(ut::complex pos, ut::complex neg, double) {
            if( std::abs(pos - neg) > 1.e-12*(std::abs(pos) + std::abs(neg)) )
                throw std::runtime_error("bath::Fit<double>: hybridisation function is not real !");
        }
    };
    
    
    template<typename HybVal>
    struct Simple {
        Simple() = delete;
        Simple(jsx::value const& jParams, std::vector<std::complex<double>> const& hyb, std::vector<std::complex<double>> hybTransp) :
        I_(std::max(static_cast<int>((jParams.is("hybridisation factor") ? jParams("hybridisation factor").real64() : 4.)*hyb.size()), 1)),
        fact_(I_/static_cast<double>(ut::KeyMax)),
        data_(I_ + 2) {
            Fit<HybVal> fit(ut::beta(), hyb, hybTransp);
            
            ut::complex const fact = -fit.moment()/(1. + std::exp((fit.eps().real() < .0 ? 1. : -1.)*fit.eps()*ut::beta()));
            for(std::size_t i = 0; i < I_ + 1; ++i) {
                double const time = ut::beta()*i/static_cast<double>(I_);
                data_[i] = to_value(fact*std::exp(fit.eps()*(fit.eps().real() < .0 ? ut::beta() - time : -time)), HybVal());
            }
            
            for(std::size_t m = 0; m < hyb.size(); ++m) {
                ut::complex const iw{.0, M_PI*(2*m + 1)/ut::beta()};
                double const smooth = (m + fit.N() < hyb.size() ? 1. : (hyb.size() - m)/static_cast<double>(fit.N()));
                ut::complex const value = (hyb[m] - fit.moment() /(iw - fit.eps()) )*smooth;
                ut::complex const valueTransp = (hybTransp[m] - std::conj(fit.moment())/(iw - std::conj(fit.eps())))*smooth;
                
                for(std::size_t i = 0; i < I_ + 1; ++i) {
                    double const time = ut::beta()*i/static_cast<double>(I_);
                    ut::complex exp{std::cos(-time*iw.imag()), std::sin(-time*iw.imag())};
                    data_[i] += to_value((exp*value + std::conj(exp*valueTransp))/ut::beta(), HybVal());
                }
            }
            
            data_.back() = .0;
        };
        Simple(Simple const&) = delete;
        Simple(Simple&&) = default;
        Simple& operator=(Simple const&) = delete;
        Simple& operator=(Simple&&) = delete;
        
        HybVal get(ut::KeyType key) const {
            double it = fact_*key; int i0 = static_cast<int>(it);
            return (1. - (it - i0))*data_[i0] + (it - i0)*data_[i0 + 1];
        };
        
        ~Simple() = default;
    private:
        std::size_t const I_;
        double const fact_;
        
        std::vector<HybVal> data_;
        
        double             to_value(ut::complex const& arg, double) const { return arg.real();};
        ut::complex const& to_value(ut::complex const& arg, ut::complex) const { return arg;};
    };
    
    
    template<typename HybVal>
    struct Hyb : itf::Hyb {
        Hyb() = delete;
        Hyb(jsx::value const& jParams, jsx::value& jMatrix, jsx::value& jFunctions) :
        flavors_(2*jParams("hybridisation")("matrix").size()),  //!!!!!!!!!! test against jAtomic !!!
        matrix_(flavors_*flavors_, nullptr) {
            mpi::cout << "Reading in hybridisation ... " << std::flush;

            for(auto const& row : jMatrix.array())
                if(row.size() != jMatrix.size())
                    throw std::runtime_error("Hyb: hybridisation is not a matrix.");

            ga::Join labels(jMatrix.size());
            for(std::size_t i = 0; i < jMatrix.size(); ++i)
                for(std::size_t j = 0; j < jMatrix.size(); ++j)
                    if(jMatrix(i)(j).string() != "") {
                        if(jMatrix(j)(i).string() == "")
                            throw std::runtime_error("Hyb: invalid hybridisation matrix.");
                        
                        std::string const entry = jMatrix(i)(j).string();
                        std::string const entryTransp = jMatrix(j)(i).string();

                        if(!data_.count(entry))
                            data_.emplace(entry, Simple<HybVal>(jParams, jsx::at<io::cvec>(jFunctions(entry)), jsx::at<io::cvec>(jFunctions(entryTransp))));

                        matrix_[(2*i + 1)  + flavors_*2*j] = &data_.at(entry);
                        
                        labels.join(i, j);
                    }
            blocks_.resize(labels.clean());

            for(std::size_t i = 0; i < jMatrix.size(); ++i)
                for(std::size_t j = 0; j < jMatrix.size(); ++j) {
                    if(labels.label(i) == labels.label(j) && jMatrix(i)(j).string() == "")
                        throw std::runtime_error("Blocks: hybridisation matrix is not a block-diagonal matrix.");
                    if(labels.label(i) != labels.label(j) && jMatrix(i)(j).string() != "")
                        throw std::runtime_error("Blocks: hybridisation matrix is not a block-diagonal matrix.");
                }

            for(std::size_t i = 0; i < jMatrix.size(); ++i) {
                blocks_[labels.label(i)].flavorsL().push_back(2*i + 1);
                blocks_[labels.label(i)].flavorsR().push_back(2*i    );
            }

            mpi::cout << "Ok" << std::endl;
        };
        Hyb(Hyb const&) = delete;
        Hyb(Hyb&&) = delete;
        Hyb& operator=(Hyb const& hyb) = delete;
        Hyb& operator=(Hyb&&) = delete;
        ~Hyb() = default;
        
        std::vector<Block> const& blocks() const { return blocks_;};
        
        HybVal operator()(int flavorL, int flavorR, ut::KeyType key) const {
            return key < 0 ? -matrix_[flavorL + flavors_*flavorR]->get(key + ut::KeyMax) : matrix_[flavorL + flavors_*flavorR]->get(key);
        };
        
    private:
        int const flavors_;
        std::map<std::string, Simple<HybVal>> data_;
        std::vector<Simple<HybVal> const*> matrix_;

        std::vector<Block> blocks_;
    };
    
    
    template<typename HybVal> Hyb<HybVal>& get(itf::Hyb& hybItf) {
        return static_cast<Hyb<HybVal>&>(hybItf);
    };
    
    template<typename HybVal> Hyb<HybVal> const& get(itf::Hyb const& hybItf) {
        return static_cast<Hyb<HybVal> const&>(hybItf);
    };
}

#endif
