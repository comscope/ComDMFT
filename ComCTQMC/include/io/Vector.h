#ifndef INCLUDE_IO_VECTOR_H
#define INCLUDE_IO_VECTOR_H

#include <vector>
#include <complex>

#include "Base64.h"
#include "../JsonX.h"

//scheiss b64 member variable !! Kack loesig im moment

namespace io {
    
    template<typename T, typename std::enable_if<!std::is_same<std::complex<double>, T>::value, int>::type = 0>   // sfinae so that it is taken for ivec and rvec. However, not sure anymore why function overload below is not sufficient
    inline jsx::value encode(std::vector<T> const& source, bool b64) {
        if(b64) return base64::encode(source);
        
        return jsx::array_t(source.begin(), source.end());
    };
    
    template<typename T, typename std::enable_if<!std::is_same<std::complex<double>, T>::value, int>::type = 0>
    inline void decode(jsx::value const& source, std::vector<T>& dest) {
        if(source.is<jsx::array_t>()) {
            for(auto const& x : source.array()) dest.push_back(x.real64());
        } else if(source.is<jsx::string_t>()) {
            base64::decode(source.string(), dest);
        } else
            throw std::runtime_error("io::decode: invalid format");
    };
    
    inline jsx::value encode(std::vector<std::complex<double>> const& source, bool b64) {
        jsx::value jDest;
        std::vector<double> real, imag;
        
        for(auto const& z : source) {
            real.push_back(z.real()); imag.push_back(z.imag());
        };
        
        jDest["real"] = encode(real, b64); jDest["imag"] = encode(imag, b64);
        
        return jDest;
    };

    
    inline void decode(jsx::value const& source, std::vector<std::complex<double>>& dest) {
        std::vector<double> real; decode(source("real"), real);
        std::vector<double> imag; decode(source("imag"), imag);
        
        if(real.size() != imag.size()) throw std::runtime_error("io::decode: invalid format");
        
        for(std::size_t n = 0; n < real.size(); ++n)
            dest.push_back({real[n], imag[n]});
    };

    
    template<typename T> struct Vector : std::vector<T> {
        inline static std::string name() { return name(T());};
        
        Vector() = default;
        Vector(Vector const&) = default;
        Vector(Vector&&) = default;
        Vector& operator=(Vector const&) = default;
        Vector& operator=(Vector&&) = default;
        ~Vector() = default;
        
        template<typename Arg, typename... Args, typename std::enable_if< !std::is_same<typename std::decay<T>::type, Vector>::value, int>::type = 0>
        Vector(Arg&& arg, Args&&... args) : std::vector<T>(std::forward<Arg>(arg), std::forward<Args>(args)...) {}
        
        template<typename Arg, typename... Args>
        Vector(std::initializer_list<Arg>&& list, Args&&... args) : std::vector<T>(std::forward<std::initializer_list<Arg>>(list), std::forward<Args>(args)...) {}
        
        template<typename Arg, typename std::enable_if< !std::is_same<typename std::decay<Arg>::type, Vector>::value, int>::type = 0>
        Vector& operator=(Arg&& arg) { std::vector<T>::operator=(std::forward<Arg>(arg)); return *this;}
        
        bool& b64() const { return b64_;};
        void read(jsx::value const& source) { decode(source, *this);};
        void write(jsx::value& dest) const { dest = encode(*this, b64_);};
    private:
        mutable bool b64_ = false;
        
        inline static std::string name(jsx::int64_t const&) { return "io::ivec";};
        inline static std::string name(double const&) { return "io::rvec";};
        inline static std::string name(std::complex<double> const&) { return "io::cvec";};
    };
    

    typedef Vector<jsx::int64_t> ivec;
    typedef Vector<double> rvec;
    typedef Vector<std::complex<double>> cvec;

};

#endif
