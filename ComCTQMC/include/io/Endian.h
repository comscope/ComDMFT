#ifndef INCLUDE_IO_ENDIAN_H
#define INCLUDE_IO_ENDIAN_H

#include <utility>
#include <cstdlib>
#include <stdexcept>
#include <vector>

namespace endian {
    
    template<typename T, bool isIntegral> struct KeyImpl {
    };
    
    template<typename T> struct KeyImpl<T, true> {
        static T get(std::size_t size = 0) { return size + 1 < sizeof(T) ? (get(size + 1) << 8) + size : size;}; // 0x....03020100 in (logical) hex
    };
    
    template<> struct KeyImpl<double, false> {
        static double get() { return 7.9499288951273625362e-275;}; // 0x0706050403020100 in (logical  sign-exponent-mantissa) hex
    };
    
    template<> struct KeyImpl<float, false> {
        static float get()  { return 3.8204714345e-37;};  // 0x03020100 in (logical  sign-exponent-mantissa) hex
    };
    
    template<typename T> struct Key{
        static T get() { return KeyImpl<T, std::is_integral<T>::value>::get();};
    };
    
    
    template<typename T>
    struct Little {
        Little() {
            std::vector<int> test(sizeof(T), 0);
            for(std::size_t index = 0; index < sizeof(T); ++index) {
                if(!( map(index) < sizeof(T) ) || test[map(index)]) throw std::runtime_error("endian::Little: something is wrong ...");
                test[map(index)] = 1;
            }
        };
        
        void write(T& arg) const {
            T temp = arg;
            
            char* arg_ptr = reinterpret_cast<char*>(&arg);
            char* tmp_ptr = reinterpret_cast<char*>(&temp);
            
            for(std::size_t i = 0; i < sizeof(T); ++i) arg_ptr[map(i)] = tmp_ptr[i];
        };
        
        void read(T& arg) const {
            T temp = arg;
            
            char* arg_ptr = reinterpret_cast<char*>(&arg);
            char* tmp_ptr = reinterpret_cast<char*>(&temp);
            
            for(std::size_t i = 0; i < sizeof(T); ++i) arg_ptr[i] = tmp_ptr[map(i)];
        };
        
    private:
        T const key_ = Key<T>::get();
        
        std::size_t map(std::size_t index) const {
            return reinterpret_cast<unsigned char const*>(&key_)[index];
        };
    };
};

#endif
