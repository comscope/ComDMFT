#ifndef UPDATES
#define UPDATES

#include <vector>
#include <complex>
#include <fstream>
#include <valarray>
#include <cstring>
#include <random>

#include "Utilities.h"

//Achtung: es kann sein dass gewisse observabeln nicht gespeichert wurden, c.f. MonteCarlo.h

namespace up {
    
    template<template<typename> class Upd, typename... Us> struct Tuple {};
    
    template<template<typename> class Upd, typename U, typename... Us> struct Tuple<Upd, U, Us...> {
        Tuple<Upd, Us...> next_;
        Upd<U> upd_;
    };
    
    template<template<typename> class Upd, typename T, typename U, typename... Us>
    struct get_impl {
        static Upd<T>& at(Tuple<Upd, U, Us...>& t) { return get_impl<Upd, T, Us...>::at(t.next_);};
    };
    
    template<template<typename> class Upd, typename T, typename... Us>
    struct get_impl<Upd, T, T, Us...> {
        static Upd<T>& at(Tuple<Upd, T, Us...>& t) { return t.upd_;};
    };
    
    template<typename T, template<typename> class Upd, typename... Us>
    Upd<T>& get(Tuple<Upd, Us...>& t) {
        return get_impl<Upd, T, Us...>::at(t);
    };
    
    //--------------------------------------------------------------------------------------------------------------------------

    struct InsertTwo {
        InsertTwo(int flavorL, int flavorR, int bath) :
        flavorL(flavorL), flavorR(flavorR), bath(bath) {
        };
        
        int const flavorL;
        int const flavorR;
        int const bath;
        
        ut::KeyType& keyL() { return  keyL_;};
        ut::KeyType& keyR() { return  keyR_;};
        
        ut::KeyType keyL() const { return  keyL_;};
        ut::KeyType keyR() const { return  keyR_;};
        
    private:
        ut::KeyType keyL_, keyR_;
    };
    
    //--------------------------------------------------------------------------------------------------------------------------
    
    struct EraseTwo {
        EraseTwo(int flavorL, int flavorR, int bath) :
        flavorL(flavorL), flavorR(flavorR), bath(bath) {
        };
        
        int const flavorL;
        int const flavorR;
        int const bath;
        
        ut::KeyType& keyL() { return  keyL_;};
        ut::KeyType& keyR() { return  keyR_;};
        
        ut::KeyType keyL() const { return  keyL_;};
        ut::KeyType keyR() const { return  keyR_;};
    
    private:
        ut::KeyType keyL_, keyR_;
    };
    
    //--------------------------------------------------------------------------------------------------------------------------
}


#endif
