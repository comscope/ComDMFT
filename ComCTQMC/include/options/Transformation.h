#ifndef INCLUDE_OPTIONS_TRANSFORMATION_H
#define INCLUDE_OPTIONS_TRANSFORMATION_H

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <bitset>
#include <cassert>
#include <iomanip>
#include <set>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <string>

#include "../JsonX.h"


namespace opt {
    
    template<typename Value>
    struct Transformation {
        Transformation(int const N, jsx::value jTransformation) :
        t_(jTransformation.is<jsx::empty_t>() ? io::PrettyMatrix<Value>(N, N) : jsx::at<io::PrettyMatrix<Value>>(jTransformation)) {
            if(jTransformation.is<jsx::empty_t>())
                for(int f = 0; f < N; ++f) t_(f, f) = 1.;
        }
        
        Value operator()(int i, int j) const {
            return t_(i, j);
        };
        
        int I() const { return t_.I();};
        int J() const { return t_.J();};
        
    private:
        io::PrettyMatrix<Value> t_;
    };
    
};

#endif






