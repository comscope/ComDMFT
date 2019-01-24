#ifndef IO
#define IO

#include <iostream>
#include <cstring>
#include <vector>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cerrno>
#include <utility>
#include <valarray>
#include <limits>
#include <map>
#include <complex>

#include "BlasLapack.h"
#include "basen.hpp"
#include "JsonX.h"

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// little-big endian madness !!!!!!! Little fish big fish swimming in the water ....
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


//scheiss b64 member variable !! Kack loesig im moment ...

namespace io {
    
    inline void encode(std::vector<double> const& source, jsx::value& dest) {
        std::string encoded;
        bn::encode_b64(reinterpret_cast<char const*>(source.data()), reinterpret_cast<char const*>(source.data() + source.size()), back_inserter(encoded));
        dest = std::move(encoded);
    }
    
    inline void decode(jsx::value const& source, std::vector<double>& dest) {
        std::string buffer; bn::decode_b64(source.string().begin(), source.string().end(), back_inserter(buffer));
        if(buffer.size()%sizeof(double) != 0) throw std::runtime_error("IO::decode: invalid data size");
        dest.resize(buffer.size()/sizeof(double)); buffer.copy(reinterpret_cast<char*>(dest.data()), buffer.size());
    }

    template<typename T>
    struct Vector : std::vector<T> {
        Vector() = default;
        Vector(Vector const&) = default;
        Vector(Vector&&) = default;
        Vector& operator=(Vector const&) = default;
        Vector& operator=(Vector&&) = default;
        
        bool& b64() const { return b64_;};
        inline void read(jsx::value const& arg) {
            throw std::runtime_error("IO::Vector::read: not implemented");
        };
        inline void write(jsx::value& arg) const {
            throw std::runtime_error("IO::Vector::write: not implemented");
        };
    private:
        mutable bool b64_ = false;
    };
    
    template<>
    inline void Vector<double>::write(jsx::value& arg) const {
        if(b64_)
            encode(*this, arg);
        else
            arg = jsx::value(this->begin(), this->end());
    };
    
    template<>
    inline void Vector<double>::read(jsx::value const& arg) {
        if(arg.type() == jsx::type::array) {
            for(auto const& x : arg.array()) this->push_back(x.real64());
        } else if(arg.type() == jsx::type::string) {
            decode(arg, *this);
        } else
            throw std::runtime_error("IO::Vector<double>::read: invalid format");
    };
    
    template<>
    inline void Vector<std::complex<double>>::write(jsx::value& arg) const {
        std::vector<double> real, imag;
        for(auto const& x : *this) {
            real.push_back(x.real());
            imag.push_back(x.imag());
        }
        if(b64_) {
            encode(real, arg["real"]);
            encode(imag, arg["imag"]);
        } else {
            arg["real"] = jsx::value(real.begin(), real.end());
            arg["imag"] = jsx::value(imag.begin(), imag.end());
        }
    };
    
    template<>
    inline void Vector<std::complex<double>>::read(jsx::value const& arg) {
        if(arg("real").type() == jsx::type::array && arg("imag").type() == jsx::type::array) {
            auto const& real = arg("real").array();
            auto const& imag = arg("imag").array();
            if(real.size() != imag.size())
                throw std::runtime_error("IO::Vector<std::complex<double>>::read: invalid format");
            for(std::size_t n = 0; n < real.size(); ++n)
                this->push_back(std::complex<double>(real[n].real64(), imag[n].real64()));
        } else if(arg("real").type() == jsx::type::string && arg("imag").type() == jsx::type::string) {
            std::vector<double> real; decode(arg("real"), real);
            std::vector<double> imag; decode(arg("imag"), imag);
            if(real.size() != imag.size())
                throw std::runtime_error("IO::Vector<std::complex<double>>::read: invalid format");
            for(std::size_t n = 0; n < real.size(); ++n)
                this->push_back(std::complex<double>(real[n], imag[n]));
        } else
            throw std::runtime_error("IO::Vector<std::complex<double>>::read: invalid format");
    };

    typedef Vector<double> rvec;
    typedef Vector<std::complex<double>> cvec;
    
    
    template<typename T>
    struct Matrix : protected Vector<T> {
        Matrix() = default;
        Matrix(Matrix const&) = default;
        Matrix(Matrix&& arg) noexcept : Vector<T>(std::move(arg)), I_(arg.I_), J_(arg.J_) { arg.I_ = arg.J_ = 0;};
        Matrix& operator=(Matrix const&) = default;
        Matrix& operator=(Matrix&& arg) { Vector<T>::operator=(std::move(arg)); I_ = arg.I_; J_ = arg.J_; arg.I_ = arg.J_ = 0; return *this;};
        
        bool& b64() const { return Vector<T>::b64();};
        
        int const& I() const { return I_;};
        int const& J() const { return J_;};
        
        T* data() { return std::vector<T>::data();};
        T const* data() const { return std::vector<T>::data();};
        
        T& operator()(int i, int j) { return this->operator[](i + j*I_);};
        T const& operator()(int i, int j) const { return this->operator[](i + j*I_);};
        
        Matrix& resize(int I, int J) {
            I_ = I; J_ = J; std::vector<T>::resize(I_*J_, .0);
            return *this;
        };
        Matrix& transpose() {
            Matrix temp; temp. std::vector<T>::resize(I_*J_);
            temp.I_ = J_; temp.J_ = I_;
            
            for(int i = 0; i < I_; ++i)
                for(int j = 0; j < J_; ++j)
                    temp(j, i) = this->operator()(i, j);
            
            return *this = std::move(temp);
        };
        
        inline void read(jsx::value const& arg) {
            throw std::runtime_error("IO::Matrix::read: not implemented");
        };
        inline void write(jsx::value& arg) const {
            throw std::runtime_error("IO::Matrix::write: not implemented");
        };
    private:
        int I_ = 0;
        int J_ = 0;
    };
    
    template<>
    inline void Matrix<double>::read(jsx::value const& arg) {
        I_ = arg(0).int64();
        J_ = arg(1).int64();
        Vector<double>::read(arg(2));
    }
    
    template<>
    inline void Matrix<double>::write(jsx::value& arg) const {
        arg = jsx::array(3);
        arg(0) = static_cast<std::int64_t>(I_);
        arg(1) = static_cast<std::int64_t>(J_);
        Vector<double>::write(arg(2));
    }
    
    typedef Matrix<double> rmat;
    typedef Matrix<std::complex<double>> cmat;
}

#endif 
