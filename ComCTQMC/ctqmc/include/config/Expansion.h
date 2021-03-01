#ifndef CTQMC_INCLUDE_CONFIG_EXPANSION_H
#define CTQMC_INCLUDE_CONFIG_EXPANSION_H

#include <stdexcept>
#include <iostream>
#include <vector>

#include "../Utilities.h"
#include "../../../include/mpi/Utilities.h"


namespace cfg {

    
    struct Entries : std::vector<ut::KeyType> {
        Entries() = default;
        Entries(Entries const&) = delete;
        Entries(Entries&&) = default; // ?????
        Entries& operator=(Entries const&) = delete;
        Entries& operator=(Entries&&) = default;
        ~Entries() = default;
        
        void insert(ut::KeyType key) {
            std::vector<ut::KeyType>::insert(std::upper_bound(begin(), end(), key), key);
        };
        
        void erase(ut::KeyType key) {
            std::vector<ut::KeyType>::erase(std::lower_bound(begin(), end(), key)); 
        };
    };
    
    
    struct Expansion : std::vector<Entries> {
        Expansion() = delete;
        Expansion(jsx::value jExpansion, std::size_t flavors) {
            if(!jExpansion.is<jsx::null_t>()) {
                if(flavors != jExpansion.size())
                    throw std::runtime_error("Expansion: invalid number of entries.");
                
                for(auto& jEntries : jExpansion.array()) {
                    Entries entries;
                    for(auto& key : jsx::at<io::ivec>(jEntries)) entries.push_back(key);
                    push_back(std::move(entries));
                }
            } else
                resize(flavors);
        };
        Expansion(Expansion const&) = delete;
        Expansion(Expansion&&) = delete;
        Expansion& operator=(Expansion const&) = delete;
        Expansion& operator=(Expansion&&) = delete;
        ~Expansion() = default;
        
        jsx::value json() const {
            jsx::array_t jExpansion;
            
            for(auto const& entries : *this) {
                io::ivec keys; keys.b64() = true;
                for(auto const& entry : entries) keys.push_back(entry);
                jExpansion.push_back(std::move(keys));
            }
            
            return jExpansion;
        };
    };

}

#endif
