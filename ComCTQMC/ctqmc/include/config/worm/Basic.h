#ifndef CTQMC_INCLUDE_CONFIG_WORM_BASIC_H
#define CTQMC_INCLUDE_CONFIG_WORM_BASIC_H


#include "../../Utilities.h"
#include "../../../../include/JsonX.h"


namespace data {
    
    template<typename Value> struct Data;
    
}

namespace state {
    
    template<typename Value> struct State;
    
}


namespace cfg {
    
    
    struct FermionicTime {
        ut::KeyType key() const { return key_;};
        ut::KeyType& key() { return key_;};
    protected:
        ut::KeyType key_;
        
        FermionicTime() = default;
        FermionicTime(ut::KeyType key) : key_(key) {};
    };
    
    struct BosonicTime {
        ut::KeyType key() const { return key_;};
        ut::KeyType& key() { return key_;};
    protected:
        ut::KeyType key_;
        
        BosonicTime() = default;
        BosonicTime(ut::KeyType key) : key_(key) {};
    };
    
    struct Flavor {
        Flavor() = default;
        Flavor(int flavor) : flavor_(flavor) {};
        
        int flavor() const { return flavor_;};
        int& flavor() { return flavor_;};
    private:
        int flavor_;
    };

    
    bool operator<(Flavor const& lhs, Flavor const& rhs) {   // this guy is dangerous ....
        return lhs.flavor() < rhs.flavor();
    };
    
    
}

#endif

