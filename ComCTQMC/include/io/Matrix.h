#ifndef IOMATRIX
#define IOMATRIX

#include <vector>
#include <complex>

#include "Vector.h"
#include "../JsonX.h"

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// little-big endian madness !!!!!!! Little fish big fish swimming in the water ....
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
        Matrix& transpose() {
            Matrix temp; temp.data_.resize(I_*J_);
            temp.I_ = J_; temp.J_ = I_;
            
            for(int i = 0; i < I_; ++i)
                for(int j = 0; j < J_; ++j)
                    temp(j, i) = operator()(i, j);
            
            return *this = std::move(temp);
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
};

#endif 
