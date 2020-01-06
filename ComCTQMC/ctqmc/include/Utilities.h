#ifndef UTILITIES_H
#define UTILITIES_H

#include <vector>
#include <iostream>
#include <complex>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <random>
#include <limits>
#include <stdexcept>


namespace ut {

	typedef std::complex<double> complex;

	typedef std::int64_t KeyType;
	KeyType const KeyMax = (static_cast<KeyType>(1) << 62);	
	inline KeyType cyclic(KeyType key) { return key < 0 ? key + KeyMax : key;}
	
    template<typename E, typename D>
    struct RandomNumberGenerator {
        RandomNumberGenerator(E const& eng, D const& distr) : eng_(eng), distr_(distr) {};
        typename D::result_type operator()() { return distr_(eng_);};
    private:
        E eng_; D distr_;
    };
    
    typedef std::mt19937 Engine;
    typedef std::uniform_real_distribution<double> UniformDistribution;
    typedef RandomNumberGenerator<Engine, UniformDistribution> UniformRng;
    
    //------------------------------ dae kack choert nit dahere ... ------------------------------------------------------------
    
    struct out_of_memory {};
    
    //--------------------------------------------------------------------------------------------------------------------------
    
    // fruender oder spoeter denn mal beta == 1, aber bis denn halt dae haesslich (aber sicheri) scheiss ....
    struct Beta {
        Beta() : beta_(nullptr) {};
        Beta(double beta) : beta_(new double) { *beta_ = beta;};
        Beta(Beta const&) = delete;
        Beta(Beta&&) = delete;
        Beta& operator=(Beta const&) = delete;
        Beta& operator=(Beta&& other) {
            beta_ = other.beta_; other.beta_ = nullptr; return *this;
        };
        ~Beta() {
            delete beta_;
        };
        
        double operator()() const { return *beta_;};

    private:
        double* beta_;
    };
    
    extern Beta beta;
    
    //--------------------------------------------------------------------------------------------------------------------------
    
    enum class Flag { Pending, Accept, Reject };
    
    //--------------------------------------------------------------------------------------------------------------------------
    
    struct Zahl {
        Zahl() : mantissa_(.0), exponent_(std::numeric_limits<int>::min()) {};
        Zahl(double x, double y = .0) {   // constructs x*exp(y)
            if(std::isfinite(x) && std::isfinite(y)) {
                mantissa_ = std::frexp(x*std::exp(y - M_LN2*(static_cast<int>(y/M_LN2) + 1)), &exponent_);
                mantissa_ ? exponent_ += static_cast<int>(y/M_LN2) + 1 : exponent_ = std::numeric_limits<int>::min();
            } else
                throw(std::runtime_error("ut::Zahl: argument of constructor is not a number"));
        };
        Zahl(Zahl const&) = default;
        Zahl(Zahl&&) = default;
        Zahl& operator=(Zahl const&) = default;
        Zahl& operator=(Zahl&&) = default;
        ~Zahl() = default;
        
        Zahl& operator+=(Zahl const& arg) {
            int exp;
            if(exponent_ > arg.exponent_) {
                mantissa_ = std::frexp(mantissa_ + std::ldexp(arg.mantissa_, arg.exponent_ - exponent_), &exp);
                mantissa_ ? exponent_ += exp : exponent_ = std::numeric_limits<int>::min();
            } else {
                mantissa_ = std::frexp(arg.mantissa_ + std::ldexp(mantissa_, exponent_ - arg.exponent_), &exp);
                mantissa_ ? exponent_ = arg.exponent_ + exp : exponent_ = std::numeric_limits<int>::min();
            }
            return *this;
        };
        Zahl& operator*=(Zahl const& arg) {
            int exp; mantissa_ = std::frexp(mantissa_*arg.mantissa_, &exp);
            mantissa_ ? exponent_ += (arg.exponent_ + exp) : exponent_ = std::numeric_limits<int>::min();
            return *this;
        };
        Zahl& operator/=(Zahl const& arg) {
            int exp; mantissa_ = std::frexp(mantissa_/arg.mantissa_, &exp);
            mantissa_ ? exponent_ += (-arg.exponent_ + exp) : exponent_ = std::numeric_limits<int>::min();
            return *this;
        };
        double to_double() const { return std::ldexp(mantissa_, exponent_);};
        Zahl abs() { mantissa_ = std::abs(mantissa_); return *this;};
        
        double mantissa() const { return mantissa_;};
        int exponent() const { return exponent_;};
    private:
        double mantissa_;
        int exponent_;
        
        friend int operator==(Zahl const&, Zahl const&);
        friend int operator<=(Zahl const&, Zahl const&);
        friend Zahl abs(Zahl const&);
    };
    
    inline int operator==(Zahl const& x, Zahl const& y) {
        return (x.mantissa_ == y.mantissa_)&&(x.exponent_ == y.exponent_);
    }
    inline int operator<=(Zahl const& x, Zahl const& y) {
        if(x.mantissa_*y.mantissa_ <= .0 || x.exponent_ == y.exponent_) return x.mantissa_ <= y.mantissa_;
        return x.exponent_ < y.exponent_ ? x.mantissa_ > .0 : x.mantissa_ < .0;
    }
    inline Zahl abs(Zahl const& arg) {
        Zahl temp(arg); return temp.abs();
    }
    
    inline Zahl exp(double arg) {
        return arg != -std::numeric_limits<double>::infinity() ? Zahl(1., arg) : Zahl();
    }	
    inline Zahl operator+(Zahl const& x, Zahl const& y) { 
        Zahl temp(x); temp += y; return temp;
    }
    inline Zahl operator*(Zahl const& x, Zahl const& y) { 
        Zahl temp(x); temp *= y; return temp;
    }
    inline Zahl operator/(Zahl const& x, Zahl const& y) {
        Zahl temp(x); temp /= y; return temp;
    }
    
    template<typename... Args> struct Options {};
}


#endif
