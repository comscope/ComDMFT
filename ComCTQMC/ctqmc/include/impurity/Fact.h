#ifndef CTQMC_INCLUDE_IMPURITY_FACT_H
#define CTQMC_INCLUDE_IMPURITY_FACT_H

#include <cmath>
#include <iostream>
#include <vector>

#include "../Utilities.h"


namespace imp {

    template<typename Value>
    struct Fact {
        Fact() = default;
        Fact(Fact const&) = delete;
        Fact(Fact&&) = delete;
        Fact& operator=(Fact const&) = delete;
        ~Fact() = default;

        void operator*=(Value arg) {
            factTry_ *= arg;
        };
        
        void operator/=(Value arg) {
            factTry_ /= arg;
        };
        
        double ratio() const {
            return std::abs(factTry_/fact_);
        };
        
        void accept() {
            fact_ = factTry_;
        };
        
        void reject() {
            factTry_ = fact_;
        };
        
        Value sign() const {
            return fact_/std::abs(fact_);
        };
        
    private:
        Value fact_{1.};
        Value factTry_{1.};
        
    };
    
}


#endif
