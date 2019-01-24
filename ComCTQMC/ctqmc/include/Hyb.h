#ifndef HYB
#define HYB

#include <cmath>
#include <iostream>
#include <sstream>
#include "Utilities.h"

namespace hy {

    template<class HybFunc>
    struct Hyb {
        Hyb() = delete;
        Hyb(jsx::value const& jParams, jsx::value& jMatrix, jsx::value& jFunctions) :
        flavors_(2*jParams("hybridisation")("matrix").size()),  //!!!!!!!!!! test against jAtomic !!!
        matrix_(flavors_*flavors_, nullptr) {
            std::cout << "Reading in hybridisation ... " << std::flush;

            for(auto const& row : jMatrix.array())
                if(row.size() != jMatrix.size())
                    throw std::runtime_error("Hyb: hybridisation is not a matrix.");

            block_.resize(jMatrix.size()); std::iota(block_.begin(), block_.end(), 0);
            for(std::size_t i = 0; i < jMatrix.size(); ++i) {
                for(std::size_t j = 0; j < jMatrix.size(); ++j)
                    if(jMatrix(i)(j).string() != "") {
                        if(jMatrix(i)(j).string() != jMatrix(j)(i).string())
                            throw std::runtime_error("Hyb: hybridisation matrix is not symmetric.");
                        
                        std::string const entry = jMatrix(i)(j).string();
                        
                        if(!data_.count(entry))
                            data_.emplace(entry, HybFunc(entry, jParams, jsx::at<io::cvec>(jFunctions(entry))));
                        
                        matrix_[(2*i + 1)  + flavors_*2*j] = &data_.at(entry);
                    }
                
                int min = jMatrix.size();
                for(std::size_t j = 0; j < jMatrix.size(); ++j)
                    if(jMatrix(i)(j).string() != "") min = std::min(block_[j], min);
                for(std::size_t j = 0; j < jMatrix.size(); ++j)
                    if(jMatrix(i)(j).string() != "") block_[j] = min;
            }
            
            for(std::size_t i = 0; i < jMatrix.size(); ++i) {
                for(std::size_t j = 0; j < jMatrix.size(); ++j) {
                    if(block_[i] == block_[j] && jMatrix(i)(j).string() == "") throw std::runtime_error("Blocks: hybridisation matrix is not a block-diagonal matrix.");
                    if(block_[i] != block_[j] && jMatrix(i)(j).string() != "") throw std::runtime_error("Blocks: hybridisation matrix is not a block-diagonal matrix.");
                }
            }
            
            std::set<int> distinct(block_.begin(), block_.end());
            
            blocks_ = distinct.size();
            
            int new_index = 0;
            for(auto const& old_index : distinct) {
                for(auto& val : block_)
                    if(val == old_index) val = new_index;
                ++new_index;
            }
            
            std::cout << "Ok" << std::endl;
        };
        Hyb(Hyb const&) = delete;
        Hyb(Hyb&&) = delete;
        Hyb& operator=(Hyb const& hyb) = delete;
        Hyb& operator=(Hyb&&) = delete;
        
        int blocks() const { return blocks_;};
        int block(int flavor) const { return block_.at(flavor);};
        
        double operator()(int flavorL, int flavorR, ut::KeyType key) const {
            return key < 0 ? -matrix_[flavorL + flavors_*flavorR]->get(key + ut::KeyMax) : matrix_[flavorL + flavors_*flavorR]->get(key);
        };
        
        ~Hyb() = default;
    private:
        int const flavors_;
        std::map<std::string, HybFunc> data_;
        std::vector<HybFunc const*> matrix_;
        
        int blocks_;
        std::vector<int> block_;
    };
    
    
	struct Simple {
        Simple() = delete;
        Simple(std::string name, jsx::value const& jParams, std::vector<std::complex<double>> const& hyb) :
        I_(std::max(static_cast<int>((jParams.is("hybridisation factor") ? jParams("hybridisation factor").real64() : 4.)*hyb.size()), 1)),
        fact_(I_/static_cast<double>(ut::KeyMax)),
        data_(I_ + 2) {
			if(!hyb.size()) throw std::runtime_error("Simple: no entries in hybridisation function " + name + "!");
			
		    int const MFit = std::max(static_cast<int>(hyb.size()*.1), 1);
			double ReD = .0, D2 = .0, ReiwD = 0.;
			for(std::size_t m = hyb.size() - MFit; m < hyb.size(); ++m) {
				ReD += hyb[m].real(); D2 += std::abs(hyb[m])*std::abs(hyb[m]); ReiwD += -M_PI*(2*m + 1)/ut::beta()*hyb[m].imag(); 
			}

			double const alpha = ReiwD/(MFit - ReD*ReD/D2);
			double const eps = -ReD/D2*alpha;
			double const fact = -alpha/(1. + std::exp(-std::abs(eps)*ut::beta()));
			for(int i = 0; i < I_ + 1; ++i) {
				double const time = eps > .0 ? ut::beta()*i/static_cast<double>(I_) : ut::beta()*(1. - i/static_cast<double>(I_));
				data_[i] = fact*std::exp(-std::abs(eps)*time);
			}

			for(unsigned int m = 0; m < hyb.size(); ++m) {
				ut::complex value = (hyb[m] - alpha/(ut::complex(.0, M_PI*(2*m + 1)/ut::beta()) - eps))*(m + MFit < hyb.size() ? 1. : (hyb.size() - m)/static_cast<double>(MFit));
				for(int i = 0; i < I_ + 1; ++i) {
					double arg = M_PI*(2*m + 1)*i/static_cast<double>(I_); 
					data_[i] += 2./ut::beta()*(value.real()*std::cos(arg) + value.imag()*std::sin(arg));
				}
			}
			data_.back() = .0;
		};
        Simple(Simple const&) = delete;
        Simple(Simple&&) = default;
        Simple& operator=(Simple const&) = delete;
        Simple& operator=(Simple&&) = delete;
		
        double get(ut::KeyType key) const {
			double it = fact_*key; int i0 = static_cast<int>(it);
			return (1. - (it - i0))*data_[i0] + (it - i0)*data_[i0 + 1];
		};
		
		~Simple() = default;
	private:
		int const I_;
        double const fact_;
		
        std::vector<double> data_;
	};
}

#endif