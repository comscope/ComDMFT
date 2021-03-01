#ifndef CTQMC_INCLUDE_IMPURITY_ALGEBRA_H
#define CTQMC_INCLUDE_IMPURITY_ALGEBRA_H

#include <cmath>
#include <iostream>
#include <vector>


namespace imp {
    
    namespace itf {
        
        template<typename Value>
        struct Batcher {
            virtual int is_ready() = 0;
            virtual void launch() = 0;
            virtual ~Batcher() = default;
        };
        
    };
    
    template<typename Mode, typename Value> struct Batcher;
    
    template<typename Mode, typename Value>
    Batcher<Mode, Value>& get(itf::Batcher<Value>& arg) {
        return static_cast<Batcher<Mode, Value>&>(arg);
    };
    
    
    template<typename Mode> struct Energies;
    
    template<typename Mode> struct Vector;
    
    template<typename Mode, typename Value> struct Matrix;
    
}


#endif
