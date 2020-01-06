#ifndef CONFIG_CONFIG_H
#define CONFIG_CONFIG_H

#include <stdexcept>
#include <iostream>
#include <vector>

#include "../Utilities.h"
#include "../../../include/mpi/Utilities.h"


namespace cfg {
    //----------------------------------------------------------------------------------------------------------------------------
    
    struct Entry {
        Entry() = default;
        explicit Entry(ut::KeyType key) : key_(key), ptr_(new double) {};
        Entry(Entry const&) = delete;
        Entry(Entry&&) = default;
        Entry& operator=(Entry const&) = delete;
        Entry& operator=(Entry&&) = default;
        ~Entry() = default;
        
        ut::KeyType key() const { return key_;};
        double* ptr() const { return ptr_.get();};
        
    private:
        ut::KeyType key_;
        std::unique_ptr<double> ptr_;
    };
    
    inline bool operator<(Entry const& lhs, Entry const& rhs) {
        return lhs.key() < rhs.key();
    }
    
    struct Entries : std::vector<Entry> {
        Entries() = default;
        Entries(Entries const&) = delete;
        Entries(Entries&&) = default; // ?????
        Entries& operator=(Entries const&) = delete;
        Entries& operator=(Entries&&) = default;
        ~Entries() = default;
        
        void insert(Entry&& entry) {
            std::vector<Entry>::insert(std::upper_bound(begin(), end(), entry), std::move(entry));
        };
        void erase(ut::KeyType key) {
            std::vector<Entry>::erase(std::lower_bound(begin(), end(), Entry(key))); // not optimal since double is constructed for nothing ... who cares.
        };
    };
    
    
    struct Config : std::vector<Entries> {
        Config() = delete;
        Config(std::int64_t mcId, std::size_t flavors) {
            mpi::cout << "Reading config file ... ";
            
            std::ifstream file("config_" + std::to_string(mcId) + ".json");
            if(file) {
                jsx::value jConfig; jsx::read(file, jConfig);
                
                if(flavors != jConfig.size())
                    throw std::runtime_error("Config: invalid number of entries.");
                
                for(auto& jEntries : jConfig.array()) {
                    Entries entries;
                    for(auto& key : jsx::at<io::ivec>(jEntries)) entries.push_back(Entry(key));
                    push_back(std::move(entries));
                }
            } else
                resize(flavors);
            
            mpi::cout << "Ok" << std::endl;
        };
        Config(Config const&) = delete;
        Config(Config&&) = delete;
        Config& operator=(Config const&) = delete;
        Config& operator=(Config&&) = delete;
        ~Config() = default;
        
        jsx::value json() const {
            jsx::array_t jConfig;
            
            for(auto const& entries : *this) {
                io::ivec keys; keys.b64() = true;
                for(auto const& entry : entries) keys.push_back(entry.key());
                jConfig.push_back(std::move(keys));
            }
            
            return jConfig;
        };
    };
    
}

#endif
