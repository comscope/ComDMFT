#ifndef CONFIG
#define CONFIG

#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include "Utilities.h"
#include "Bath.h"

// This needs to be modified: separate functionality (propose update) form the data (the configs) !!!!!!!
// Test low order !!!!! k = 0, 2 !!!!!!

namespace co {
	//----------------------------------------------------------------------------------------------------------------------------
    struct Abstract {
        virtual int propose(up::InsertTwo& u, ut::UniformRng& urng) = 0;
        virtual double ratio(up::InsertTwo const& u) = 0;
        virtual void accept(up::InsertTwo const& u) = 0;
        virtual void reject(up::InsertTwo const& u) = 0;
        
        virtual int propose(up::EraseTwo& u, ut::UniformRng& urng) = 0;
        virtual double ratio(up::EraseTwo const& u) = 0;
        virtual void accept(up::EraseTwo const& u) = 0;
        virtual void reject(up::EraseTwo const& u) = 0;
        
        virtual ~Abstract() {};
    };
    
    
    struct ConfigStd : Abstract {
        ConfigStd() = delete;
        ConfigStd(std::size_t flavors, std::vector<ba::Bath> const& baths) : ops_(2*flavors) {
            for(auto const& bath : baths) {
                for(auto const& op : bath.opsL()) ops_[op.flavor()].push_back(op.key());
                for(auto const& op : bath.opsR()) ops_[op.flavor()].push_back(op.key());
            }
        }
        ConfigStd(ConfigStd const&) = delete;
        ConfigStd(ConfigStd&&) = delete;
        ConfigStd& operator=(ConfigStd const&) = delete;
        ConfigStd& operator=(ConfigStd&&) = delete;
        
        int propose(up::InsertTwo& u, ut::UniformRng& urng) {   
            u.keyL() = urng()*ut::KeyMax;
            u.keyR() = urng()*ut::KeyMax;
            return 1;
        };
        double ratio(up::InsertTwo const& u) {
            return ut::beta()*ut::beta()/((ops_[u.flavorL].size() + 1.)*(ops_[u.flavorR].size() + 1.));
        };
        void accept(up::InsertTwo const& u) {
            ops_[u.flavorL].push_back(u.keyL());
            ops_[u.flavorR].push_back(u.keyR());
        };
        void reject(up::InsertTwo const& u) {};
        
        int propose(up::EraseTwo& u, ut::UniformRng& urng) {
            if(!(ops_[u.flavorL].size() && ops_[u.flavorR].size())) return 0;
            
            posL_ = urng()*ops_[u.flavorL].size(); u.keyL() = ops_[u.flavorL][posL_];
            posR_ = urng()*ops_[u.flavorR].size(); u.keyR() = ops_[u.flavorR][posR_];
            
            return 1;
        };
        double ratio(up::EraseTwo const& u) {
            return (ops_[u.flavorL].size()*ops_[u.flavorR].size())/(ut::beta()*ut::beta());
        };
        void accept(up::EraseTwo const& u) {
            ops_[u.flavorL].erase(ops_[u.flavorL].begin() + posL_);
            ops_[u.flavorR].erase(ops_[u.flavorR].begin() + posR_);
        };
        void reject(up::EraseTwo const& u) {};
    private:
        std::size_t posL_;
        std::size_t posR_;

        std::vector<std::vector<ut::KeyType>> ops_;
    };
    
    
    struct ConfigCSQ : Abstract {
        ConfigCSQ() = delete;
        ConfigCSQ(std::size_t flavors, std::vector<ba::Bath> const& baths) : ops_(2*flavors) {
            for(auto const& bath : baths) {
                for(auto const& op : bath.opsL()) insert(op.flavor(), op.key());
                for(auto const& op : bath.opsR()) insert(op.flavor(), op.key());
            }
        }
        ConfigCSQ(ConfigCSQ const&) = delete;
        ConfigCSQ(ConfigCSQ&&) = delete;
        ConfigCSQ& operator=(ConfigCSQ const&) = delete;
        ConfigCSQ& operator=(ConfigCSQ&&) = delete;
        
        int propose(up::InsertTwo& u, ut::UniformRng& urng) {
            auto const& opsL = ops_[u.flavorL];
            auto const& opsR = ops_[u.flavorR];

            ut::KeyType keyHigh;
            if(opsL.size() && opsR.size()) {
                std::size_t const index = urng()*(opsL.size() + opsR.size());
                ut::KeyType const keyLow = index < opsL.size() ? opsL[index] : opsR[index - opsL.size()];
                
                auto it = std::upper_bound(opsL.begin(), opsL.end(), keyLow);
                keyHigh = it != opsL.end() ? *it : *opsL.begin() + ut::KeyMax;
                
                it = std::upper_bound(opsR.begin(), opsR.end(), keyLow);
                keyHigh = std::min(it != opsR.end() ? *it : *opsR.begin() + ut::KeyMax, keyHigh);
                
                keyDiff_ = keyHigh - keyLow; if(keyHigh > ut::KeyMax) keyHigh -= ut::KeyMax;
            } else
                keyHigh = keyDiff_ = ut::KeyMax;
        
            u.keyL() = ut::cyclic(keyHigh - urng()*keyDiff_);
            u.keyR() = ut::cyclic(keyHigh - urng()*keyDiff_);
            
            return 1;
        };
        double ratio(up::InsertTwo const& u) {
            int const N = ops_[u.flavorL].size() + ops_[u.flavorR].size();
            double const timeDiff = (ut::beta()*keyDiff_)/ut::KeyMax;
            return timeDiff*timeDiff*(N ? N/(N + 2.) : 1.);
        };
        void accept(up::InsertTwo const& u) {
            insert(u.flavorL, u.keyL());
            insert(u.flavorR, u.keyR());
        };
        void reject(up::InsertTwo const& u) {};
        
        int propose(up::EraseTwo& u, ut::UniformRng& urng) {
            auto const& opsL = ops_[u.flavorL];
            auto const& opsR = ops_[u.flavorR];
            
            if(!(opsL.size() && opsR.size())) return 0;
            
            if(opsL.size() + opsR.size() > 2) {
                std::size_t index = urng()*(opsL.size() + opsR.size());
                ut::KeyType const keyLow = index < opsL.size() ? opsL[index] : opsR[index - opsL.size()];
                
                auto itL = std::upper_bound(opsL.begin(), opsL.end(), keyLow); ut::KeyType shiftL = 0;
                if(itL != opsL.end()) u.keyL() = *itL; else { u.keyL() = *(itL = opsL.begin()); shiftL = ut::KeyMax;};
                
                auto itR = std::upper_bound(opsR.begin(), opsR.end(), keyLow); ut::KeyType shiftR = 0;
                if(itR != opsR.end()) u.keyR() = *itR; else { u.keyR() = *(itR = opsR.begin()); shiftR = ut::KeyMax;};
                
                ut::KeyType keyHigh = ++itL != opsL.end() ? *itL + shiftL : *opsL.begin() + ut::KeyMax;
                keyHigh = std::min(++itR != opsR.end() ? *itR + shiftR : *opsR.begin() + ut::KeyMax, keyHigh);
                
                if(std::max(u.keyL() + shiftL, u.keyR() + shiftR) < keyHigh) keyDiff_ = keyHigh - keyLow; else return 0;
            } else {
                u.keyL() = *opsL.begin();
                u.keyR() = *opsR.begin();
                keyDiff_ = ut::KeyMax;
            }
            
            return 1;
        };
        double ratio(up::EraseTwo const& u) {
            int const N = ops_[u.flavorL].size() + ops_[u.flavorR].size();
            double const timeDiff = (ut::beta()*keyDiff_)/ut::KeyMax;
            return (N > 2 ? N/(N - 2.) : 1.)/(timeDiff*timeDiff);
        };
        void accept(up::EraseTwo const& u) {
            erase(u.flavorL, u.keyL());
            erase(u.flavorR, u.keyR());
        };
        void reject(up::EraseTwo const& u) {};
    private:
        ut::KeyType keyDiff_;
        std::vector<std::vector<ut::KeyType>> ops_;
        
        void insert(int flavor, ut::KeyType key) {
            auto& ops = ops_[flavor]; ops.insert(std::upper_bound(ops.begin(), ops.end(), key), key);
        };
        void erase(int flavor, ut::KeyType key) {
            auto& ops = ops_[flavor]; ops.erase(std::lower_bound(ops.begin(), ops.end(), key));
        };
    };

}

#endif