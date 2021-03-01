#ifndef INCLUDE_OPTIONS_SPHERICALHARMONICS_CLEBSCHGORDAN_H
#define INCLUDE_OPTIONS_SPHERICALHARMONICS_CLEBSCHGORDAN_H

#include <cmath>

// apply J^2 to |JM> = \sum |j1m1,j2m2><j1m1,j2m2|JM>

namespace cg {
    
    namespace impl {
        
        inline double clebschGordan(int j1, int m1, int j2, int m2, int J, int M) {
            if((j1 < 0) || (j2 < 0) || (J < 0) || ((j1 + m1)%2) || ((j2 + m2)%2) || ((J + M)%2))
                throw std::runtime_error("cg::impl::clebschGordan: invalid arguments");
            
            if((m1 + m2 != M) || (std::abs(j1 - j2) > J) || (J > j1 + j2) || ((j1 + j2 + J)%2) || (std::abs(m1) > j1) || (std::abs(m2) > j2) || (std::abs(M) > J))
                return .0;
            
            int const m1_low  = (M - j1 - j2 + std::abs(j1 - j2 + M))/2;
            int const m1_high = (M + j1 + j2 - std::abs(j1 - j2 - M))/2;
            
            double b_prev = .0, x_prev = .0, x = 1., sum = 1., coeff = 1.;
            
            for(int m = m1_high; m > m1_low; m -= 2) {
                double const a = (j1*(j1 + 2.) + j2*(j2 + 2.) + 2.*m*(M - m) - J*(J + 2.))/4.;
                double const b = std::sqrt(j1*(j1 + 2.) - (m - 2.)*m)/2.*std::sqrt(j2*(j2 + 2.) - (M - (m - 2.))*(M - m))/2.;
                double const x_temp = x;
                
                x = -(x*a + x_prev*b_prev)/b;
                sum += x*x;
                if(m - 2 == m1) coeff = x;
                
                b_prev = b; x_prev = x_temp;
            }
            
            // TODO: check if first row gives zero !
            
            return coeff/std::sqrt(sum);
        }
        
    }
    
    
    inline double clebschGordan(double j1, double m1, double j2, double m2, double J, double M) {
        return impl::clebschGordan(std::round(2.*j1), std::round(2.*m1), std::round(2.*j2), std::round(2.*m2), std::round(2.*J), std::round(2.*M));
    }
}

#endif
