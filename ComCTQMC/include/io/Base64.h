#ifndef INCLUDE_IO_BASE64_H
#define INCLUDE_IO_BASE64_H

#include <utility>
#include <cstdlib>
#include <stdexcept>
#include <vector>

#include "Endian.h"

namespace base64 {
    
    inline unsigned char get_six_at(unsigned char const* begin, std::size_t bit) {
        std::size_t const offset = bit%8;
        
        unsigned char const byte1 = begin[bit/8];
        if(offset == 0) return ((byte1     ) & 0x3F);
        if(offset == 2) return ((byte1 >> 2) & 0x3F);
        
        unsigned char const byte2 = begin[bit/8 + 1];
        if(offset == 4) return ((byte1 >> 4) & 0x0F) | ((byte2 << 4) & 0x30);
        if(offset == 6) return ((byte1 >> 6) & 0x03) | ((byte2 << 2) & 0x3C);
        
        throw std::runtime_error("base64::get_six_at: invalid position"); 
    };
    
    
    inline void set_six_at(unsigned char pattern, unsigned char* begin, std::size_t bit) {
        std::size_t const offset = bit%8;
        
        unsigned char& byte1 = begin[bit/8];
        if(offset == 0) { byte1 = (byte1 & 0xC0) | ( pattern       & 0x3F); return;};
        if(offset == 2) { byte1 = (byte1 & 0x03) | ((pattern << 2) & 0xFC); return;};
        
        unsigned char& byte2 = begin[bit/8 + 1];
        if(offset == 4) { byte1 = (byte1 & 0x0F) | ((pattern << 4) & 0xF0); byte2 = (byte2 & 0xFC) | ((pattern >> 4) & 0x03); return;};
        if(offset == 6) { byte1 = (byte1 & 0x3F) | ((pattern << 6) & 0xC0); byte2 = (byte2 & 0xF0) | ((pattern >> 2) & 0x0F); return;};
        
        throw std::runtime_error("base64::set_six_at: invalid position");
    };
    
    
    struct Dictionary {
        Dictionary() : decode_(256, 64) {
            for(unsigned char six_bits = 0; six_bits < 64; ++six_bits) {
                auto& entry = decode_[*reinterpret_cast<unsigned char const*>(&encode_[six_bits])];
                if(entry != 64) throw std::runtime_error("base64::constructor: invalid dictionary");
                entry = six_bits;
            }
        };
        ~Dictionary() = default;
        
        char encode(unsigned char arg) const {
            if(arg & 0xC0) throw std::runtime_error("base64::Dictionary::encode: invalid key");
            return encode_[arg];
        };
        
        unsigned char decode(char arg) const {
            auto temp = decode_[*reinterpret_cast<unsigned char*>(&arg)];
            if(temp == 64) throw std::runtime_error("base64::Dictionary::decode: invalid key");
            return temp;
        };
        
    private:
        char const* const encode_ = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        std::vector<unsigned char> decode_;
    };
    
    
    template<typename T>
    std::string encode(std::vector<T> const& source) {
        std::string dest;
        
        auto temp = source; endian::Little<T> as_little; for(auto& x : temp) as_little.write(x);
        
        temp.push_back(T()); auto const begin = reinterpret_cast<unsigned char const*>(temp.data()); Dictionary dict;
        
        for(std::size_t bit = 0; bit < 8*sizeof(T)*source.size(); bit += 6) dest.push_back(dict.encode(get_six_at(begin, bit)));
        
        return dest;
    };
    
    
    template<typename T>
    void decode(std::string const& source, std::vector<T>& dest) {
        dest.resize((6*source.size())/(8*sizeof(T)));
        
        if(!(6*source.size() - 8*sizeof(T)*dest.size() < 6)) throw std::runtime_error("base64::decode: too much padding bits");
        
        dest.push_back(T()); auto const begin = reinterpret_cast<unsigned char*>(dest.data()); Dictionary dict;
        
        for(std::size_t index = 0; index < source.size(); ++index) set_six_at(dict.decode(source[index]), begin, 6*index);
        
        dest.pop_back(); endian::Little<T> as_little; for(auto& x : dest) as_little.read(x);
    };
};

#endif
