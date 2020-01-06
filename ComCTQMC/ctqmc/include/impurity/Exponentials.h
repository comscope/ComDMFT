#ifndef IMPURITY_EXPONENTIALS_H
#define IMPURITY_EXPONENTIALS_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <valarray>

#include "../config/Config.h"
#include "../Utilities.h"
#include "../Data.h"

namespace imp {
    
    struct Exponentials {
        Exponentials() = delete;
        Exponentials(jsx::value const& jParams, data::Data const& data) :
        flavors_(data.ops().flavors()/2),
        nMat_(std::max(static_cast<int>(ut::beta()*(jParams.is("susceptibility cutoff") ? jParams("susceptibility cutoff").real64() : .0)/(2*M_PI)), 1)),
        ops_(2*flavors_),
        data_(flavors_, std::valarray<ut::complex>(.0, nMat_)) {
        };
        Exponentials(Exponentials const&) = delete;
        Exponentials(Exponentials&&) = delete;
        Exponentials& operator=(Exponentials const&) = delete;
        Exponentials& operator=(Exponentials&&) = delete;
        ~Exponentials() = default;
        
        void set(cfg::Config const& config) {
            std::vector<std::vector<ut::KeyType>> insert(2*flavors_), temp(2*flavors_);
            
            for(int flavor = 0; flavor < 2*flavors_; ++flavor)
                for(auto const& c : config[flavor]) {
                    temp[flavor].push_back(c.key()); auto op = std::find(ops_[flavor].begin(), ops_[flavor].end(), c.key());
                    if(op == ops_[flavor].end()) insert[flavor].push_back(c.key()); else ops_[flavor].erase(op);
                }
            
            std::vector<std::vector<ut::KeyType>> erase = std::move(ops_); ops_ = std::move(temp);
            
            for(int flavor = 0; flavor < flavors_; ++flavor)
                if(insert[2*flavor].size() + erase[2*flavor].size() + insert[2*flavor + 1].size() + erase[2*flavor + 1].size() < ops_[2*flavor].size() + ops_[2*flavor + 1].size()) {
                    update(insert[2*flavor], data_[flavor],  1); update(insert[2*flavor + 1], data_[flavor], -1);
                    update( erase[2*flavor], data_[flavor], -1); update( erase[2*flavor + 1], data_[flavor],  1);
                } else {
                    data_[flavor] = .0;
                    update(ops_[2*flavor], data_[flavor], 1); update(ops_[2*flavor + 1], data_[flavor], -1);
                }
        }
        
        void clean(cfg::Config const& config) {
            for(int flavor = 0; flavor < flavors_; ++flavor) {
                ops_[2*flavor].clear();
                for(auto const& c : config[2*flavor]) ops_[2*flavor].push_back(c.key());
                
                ops_[2*flavor + 1].clear();
                for(auto const& c : config[2*flavor + 1]) ops_[2*flavor + 1].push_back(c.key());
                
                data_[flavor] = .0;
                update(ops_[2*flavor], data_[flavor], 1); update(ops_[2*flavor + 1], data_[flavor], -1);
            }
        };
        
        std::valarray<ut::complex> const& at(int f) const { return data_[f];};
        
    private:
        int const flavors_;
        std::size_t const nMat_;
        std::vector<std::vector<ut::KeyType>> ops_;
        std::vector<std::valarray<ut::complex>> data_;
        
        void update(std::vector<ut::KeyType> const& keys, std::valarray<ut::complex>& data, int operation) {
            for(auto const& key : keys) {
                double const u = key/static_cast<double>(ut::KeyMax);
                
                ut::complex const rotate = ut::complex(std::cos(2*M_PI*u), std::sin(2*M_PI*u)); ut::complex value = operation;
                for(std::size_t n = 1; n < nMat_; ++n) {
                    value *= rotate; data[n] += value;
                }
            }
        };
    };
    
}

#endif
