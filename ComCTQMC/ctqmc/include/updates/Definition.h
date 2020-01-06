#ifndef UPDATES_DEFINITION_H
#define UPDATES_DEFINITION_H

#include <vector>
#include <cstring>
#include <cstdint>


namespace upd {
    
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
        
        double*& ptrL() { return ptrL_;};
        double*& ptrR() { return ptrR_;};
        
        double* ptrL() const { return ptrL_;};
        double* ptrR() const { return ptrR_;};
        
    private:
        ut::KeyType keyL_, keyR_;
        double *ptrL_, *ptrR_;
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
    
}


#endif
