#ifndef EVALSIM_PARTITON_BESSEL_H
#define EVALSIM_PARTITON_BESSEL_H

#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <limits>

namespace evalsim {
    
    namespace be {
        
        struct Bessel {
            Bessel() = delete;
            Bessel(std::size_t const N, std::size_t const M) : data_(N, std::vector<double>(M)) {
                for(std::size_t m = 0; m < M; ++m) {
                    std::size_t const n_max = N + 1000000; // maybe a little bit exagerated, eventually ....
                    
                    std::vector<double> result; std::vector<double> prev = bessel(N + 10, m);
                    std::size_t converged = 0;
                    
                    for(std::size_t n_high = N + 20; n_high < n_max; n_high += 10) {
                        result = bessel(n_high, m);
                        for(std::size_t n = 0; n < N; ++n)
                            if(std::abs(result[n] - prev[n]) > 1.e-14*std::abs(result[n])) converged = 0;
                        if(++converged == 2) break;
                        prev = std::move(result);
                    }
                    
                    if(converged != 2) throw std::runtime_error("Bessel: not converged !");
                    
                    for(std::size_t n = 0; n < N; ++n) data_[n][m] = result[n];
                }
            }
            Bessel(Bessel const&) = delete;
            Bessel(Bessel&&) = default;
            Bessel& operator=(Bessel const&) = delete;
            Bessel& operator=(Bessel&&) = delete;
            ~Bessel() = default;
            
            double operator()(std::size_t n, std::size_t m) {
                return data_[n][m];
            };
            
        private:
            std::vector<std::vector<double>> data_;
            
            std::vector<double> bessel_forward(std::size_t const n_high, std::size_t const m) {
                double const omega = (2.*m + 1.)*M_PI/2.;
                
                std::vector<double> result(n_high + 1);
                
                result[0] = (m%2 ? -1. : 1.)/omega; if(n_high) result[1] = 1./omega*result[0];
                for(std::size_t n = 1; n < n_high; ++n)
                    result[n + 1] = (2.*n + 1.)/omega*result[n] - result[n - 1];
                
                return result;
            };
            
            std::vector<double> bessel_backward(std::size_t const n_low, std::size_t const n_high, std::size_t const m, double const val_before_low) {
                double const omega = (2.*m + 1.)*M_PI/2.;
                double const max_val = std::numeric_limits<double>::max()/(100.*std::max((2.*n_high + 1.)/omega, 1.));
                
                std::vector<double> result(n_high - n_low + 3);
                result[n_high + 2 - n_low] = 0.; result[n_high + 1 - n_low] = 1.;
                
                for(std::size_t n = n_high + 1; n > n_low; --n) {
                    result[n - 1 - n_low] = (2.*n + 1.)/omega*result[n - n_low] - result[n + 1 - n_low];
                    if(std::abs(result[n - 1 - n_low]) > max_val) {
                        double const re_norm = 1./result[n - 1 - n_low];
                        for(std::size_t k = n - 1; k <= n_high; ++k) result[k - n_low] *= re_norm;
                    }
                }
                
                double const scal = val_before_low/((2.*n_low + 1.)/omega*result[0] - result[1]);
                for(auto& x : result) x *= scal;
                
                result.pop_back(); result.pop_back(); return result;
            };
            
            std::vector<double> bessel(std::size_t const n, std::size_t const m) {
                double const omega = (2.*m + 1.)*M_PI/2.;
                std::size_t const n_trans = static_cast<std::size_t>(omega);
                
                if(n_trans >= n) return bessel_forward(n, m);
                
                std::vector<double> forward  = bessel_forward(n_trans, m);
                std::vector<double> backward = bessel_backward(n_trans + 1, n, m, forward[n_trans]);
                
                forward.insert(forward.end(), backward.begin(), backward.end());
                
                return forward;
            };
        };
        
    }
    
}

#endif


