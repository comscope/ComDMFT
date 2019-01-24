#ifndef GREEN
#define GREEN

#include <cmath>
#include <stdexcept>
#include <iostream>
#include "Utilities.h"
#include "Bath.h"
#include "../../include/Measurements.h"

namespace gr {
    
    //viel z'viel beta's wo da umelungerend ....
    
    template<class Meas>
    struct Green {
        Green() = delete;
        Green(jsx::value const& jParams) :
        flavors_(2*jParams("hybridisation")("matrix").size()),
        matrix_(flavors_*flavors_, nullptr)
        {
            auto const& jMatrix = jParams("hybridisation")("matrix");
            for(std::size_t i = 0; i < jMatrix.size(); ++i)
                for(std::size_t j = 0; j < jMatrix.size(); ++j)
                    if(jMatrix(i)(j).string() != "") {
                        std::string const entry = jMatrix(i)(j).string();
                        
                        if(!data_.count(entry))
                            data_.emplace(entry, std::make_pair(0, Meas(jParams)));
                        
                        data_.at(entry).first += 1;
                        matrix_[2*j  + flavors_*(2*j + 1)] = &data_.at(entry).second;  //hyb isch symmetrisch, c.f. Hyb.h Muscht du aber nochmal gucken wegen allgemeinem fall ....
                    }
        };
        Green(Green const&) = delete;
        Green(Green&&) = delete;
        Green& operator=(Green const&) = delete;
        Green& operator=(Green&&) = delete;
        
        void sample(int sign, ba::Bath const& bath) {
            double const* data = bath.B().data();
            for(auto const& opL : bath.opsL())
                for(auto const& opR : bath.opsR())
                    matrix_[opR.flavor() + flavors_*opL.flavor()]->add(opR.key() - opL.key(), sign**data++);
        };
        
        void store(jsx::value& measurements, std::int64_t samples) {
            for(auto& entry : data_)
                entry.second.second.store(measurements["Green"][entry.first], entry.second.first*samples);
        };
        
        ~Green() = default;
    private:
        int const flavors_;
        std::map<std::string, std::pair<std::int64_t, Meas>> data_;
        std::vector<Meas*> matrix_;
    };
    
    
	struct BinMoments {
        BinMoments() = delete;
		BinMoments(jsx::value const& jParams) :
		nMatG_(std::max(static_cast<int>(ut::beta()*jParams("green cutoff").real64()/(2*M_PI)), 1)), 
		nItG_(4*(2*nMatG_ + 1)), 
		DeltaInv_(nItG_/ut::beta()),
		green_(4*nItG_, .0) {
        };
        BinMoments(BinMoments const&) = delete;
        BinMoments(BinMoments&&) = default;
        BinMoments& operator=(BinMoments const&) = delete;
        BinMoments& operator=(BinMoments&&) = delete;
		
        void add(ut::KeyType key, double value) {
            if(key < 0) {
                key += ut::KeyMax;
                value *= -1.;
            }
            
            double const time = key*ut::beta()/ut::KeyMax;
			int const index = static_cast<int>(DeltaInv_*time);
			double const Dtime = time - static_cast<double>(index + .5)/DeltaInv_;
			
			double* green = green_.data() + 4*index;
			*green++ += value; 
			value *= Dtime;
			*green++ += value;
			value *= Dtime;
			*green++ += value;
			value *= Dtime;
			*green += value;
		};
		
        void store(jsx::value& measurements, std::int64_t samples) {
            std::vector<std::complex<double>> green(nMatG_);
			double Dtau = ut::beta()/static_cast<double>(nItG_);
			
			for(int m = 0; m < nMatG_; ++m) {
				double omega = M_PI*static_cast<double>(2*m + 1)/ut::beta();
				double lambda = -2.*std::sin(omega*Dtau/2.)/((Dtau*omega*(1. - omega*omega*Dtau*Dtau/24.))*ut::beta());  //missing -1/beta factor
				
				ut::complex iomega(.0, omega);
				ut::complex exp(std::exp(iomega*Dtau/2.));
				ut::complex fact(std::exp(iomega*Dtau));  
				
                std::complex<double> temp = .0;
				double const* itGreen = green_.data();
				for(int i = 0; i < nItG_; ++i) {
					ut::complex coeff = lambda*exp;
					temp += coeff**itGreen++;
					coeff *= iomega;
					temp += coeff**itGreen++;
					coeff *= iomega/2.;
					temp += coeff**itGreen++;
					coeff *= iomega/3.;
					temp += coeff**itGreen++;
					
					exp *= fact;
				}
				
                green[m] = temp;
			}
			
            measurements << meas::fix(green, samples);
            
			std::memset(green_.data(), 0, 4*nItG_*sizeof(double));
		};
		~BinMoments() {};
	private:
		int const nMatG_;
		int const nItG_;
		double const DeltaInv_;
		
        std::vector<double> green_;
	};	
	
	struct OrthogonalPolynomials {
		// äs bitzeli kosmetik
	};
}
#endif