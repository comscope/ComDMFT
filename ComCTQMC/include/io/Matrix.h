#ifndef INCLUDE_IO_MATRIX_H
#define INCLUDE_IO_MATRIX_H

#include <vector>
#include <complex>

#include "Vector.h"
#include "../JsonX.h"
#include "../../ctqmc/include/Utilities.h"

//scheiss b64 member variable !! Kack loesig im moment ...

namespace io {
    
    template<typename T> struct Matrix {
        inline static std::string name() { return name(T());};
        
        Matrix() = default;
        Matrix(std::size_t I, std::size_t J, T value = .0) : I_(I), J_(J), data_(I_*J_, value) {};
        Matrix(Matrix const&) = default;
        Matrix(Matrix&& other) noexcept : I_(other.I_), J_(other.J_), data_(std::move(other.data_)) { other.I_ = other.J_ = 0;};
        Matrix& operator=(Matrix const&) = default;
        Matrix& operator=(Matrix&& other) { I_ = other.I_; J_ = other.J_; other.I_ = other.J_ = 0; data_ = std::move(other.data_);  return *this;};
        ~Matrix() = default;

        int const& I() const { return I_;};
        int const& J() const { return J_;};
        
        T* data() { return data_.data();};
        T const* data() const { return data_.data();};
        
        T& operator()(int i, int j) { return data_.at(i + j*I_);};
        T const& operator()(int i, int j) const { return data_.at(i + j*I_);};
        
        Matrix& resize(int I, int J, T value = .0) {
            I_ = I; J_ = J; data_.resize(I_*J_, value);
            return *this;
        };
        Matrix& conj() {
            Matrix temp; temp.data_.resize(I_*J_);
            temp.I_ = J_; temp.J_ = I_;
            
            for(int i = 0; i < I_; ++i)
                for(int j = 0; j < J_; ++j)
                    temp(j, i) = ut::conj(operator()(i, j));
            
            return *this = std::move(temp);
        };
        
        Matrix conj() const {
            Matrix temp; temp.data_.resize(I_*J_);
            temp.I_ = J_; temp.J_ = I_;
            
            for(int i = 0; i < I_; ++i)
                for(int j = 0; j < J_; ++j)
                    temp(j, i) = ut::conj(operator()(i, j));
            
            return temp;
        };
        
        void read(jsx::value const& source) {
            I_ = source(0).int64();
            J_ = source(1).int64();
            data_.read(source(2));
        };
        void write(jsx::value& dest) const {
            dest = jsx::array_t(3);
            
            dest(0) = static_cast<jsx::int64_t>(I());
            dest(1) = static_cast<jsx::int64_t>(J());
            data_.write(dest(2));
        };
        bool& b64() const {
            return data_.b64();
        };
        
    private:
        int I_ = 0, J_ = 0;
        Vector<T> data_;

        inline static std::string name(double const&) { return "io::rmat";};
        inline static std::string name(std::complex<double> const&) { return "io::cmat";};
    };

    
    typedef Matrix<double> rmat;
    typedef Matrix<std::complex<double>> cmat;
    
    
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     Patch    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    template<typename T> struct PrettyMatrix {
        inline static std::string name() { return name(T());};
        
        PrettyMatrix() = default;
        PrettyMatrix(std::size_t I, std::size_t J, T value = .0) : I_(I), J_(J), data_(I_*J_, value) {};
        PrettyMatrix(PrettyMatrix const&) = default;
        PrettyMatrix(PrettyMatrix&& other) noexcept : I_(other.I_), J_(other.J_), data_(std::move(other.data_)) { other.I_ = other.J_ = 0;};
        PrettyMatrix& operator=(PrettyMatrix const&) = default;
        PrettyMatrix& operator=(PrettyMatrix&& other) { I_ = other.I_; J_ = other.J_; other.I_ = other.J_ = 0; data_ = std::move(other.data_);  return *this;};
        ~PrettyMatrix() = default;
        
        int const& I() const { return I_;};
        int const& J() const { return J_;};
        
        T& operator()(int i, int j) { return data_.at(J_*i + j);};
        T const& operator()(int i, int j) const { return data_.at(J_*i + j);};
        
        PrettyMatrix& resize(int I, int J, T value = .0) {
            I_ = I; J_ = J; data_.resize(I_*J_, value);
            return *this;
        };

        void read(jsx::value const& source);
        void write(jsx::value& dest) const;
        
    private:
        int I_ = 0, J_ = 0;
        std::vector<T> data_;
        
        inline static std::string name(double const&) { return "io::prettyrmat";};
        inline static std::string name(std::complex<double> const&) { return "io::prettycmat";};
    };
    
    
    template<>
    inline void PrettyMatrix<double>::read(jsx::value const& source) {
        data_.resize(0); I_ = source.size(); J_ = I_ ? source(0).size() : 0;
        
        for(auto const& row : source.array()) {
            if(row.size() != J_) throw std::runtime_error("io::prettyrmat: wrong format");
            
            for(auto const& elem : row.array()) data_.push_back(elem.real64());
        }
    };
    
    template<>
    inline void PrettyMatrix<double>::write(jsx::value& dest) const {
        dest = jsx::array_t(I_);
        for(int i = 0; i < I_; ++i)
            dest(i) = jsx::array_t(data_.begin() + i*J_, data_.begin() + (i + 1)*J_);
    };
    
    
    template<>
    inline void PrettyMatrix<std::complex<double>>::read(jsx::value const& source) {
        if(!(source.is("real") && source.is("imag") && source.size() == 2))
            throw std::runtime_error("io::prettycmat: wrong format");
        
        data_.resize(0); I_ = source("real").size(); J_ = I_ ? source("real")(0).size() : 0;
        
        if(I_ != source("imag").size())
            throw std::runtime_error("io::prettycmat: wrong format");

        for(int i = 0; i < I_; ++i) {
            if(source("real")(i).size() != J_) throw std::runtime_error("io::PrettyMatrix: wrong format");
            if(source("imag")(i).size() != J_) throw std::runtime_error("io::PrettyMatrix: wrong format");
            
            for(int j = 0; j < J_; ++j)
                data_.push_back({source("real")(i)(j).real64(), source("imag")(i)(j).real64()});
        }
    };
    
    template<>
    inline void PrettyMatrix<std::complex<double>>::write(jsx::value& dest) const {
        jsx::value real = jsx::array_t(I_, jsx::array_t(J_));
        for(int i = 0; i < I_; ++i)
            for(int j = 0; j < J_; ++j)
                real(i)(j) = operator()(i, j).real();
        
        jsx::value imag = jsx::array_t(I_, jsx::array_t(J_));
        for(int i = 0; i < I_; ++i)
            for(int j = 0; j < J_; ++j)
                imag(i)(j) = operator()(i, j).imag();
        
        dest = jsx::object_t{{"real", real}, {"imag", imag}};
    };
    
    
    typedef PrettyMatrix<double> prettyrmat;
    typedef PrettyMatrix<std::complex<double>> prettycmat;
};

#endif 
