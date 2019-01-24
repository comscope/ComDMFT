#ifndef ZAHL
#define ZAHL

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace za {
	struct Zahl {
		Zahl(double x = .0, double y = .0) { 
			if(std::isfinite(x) && std::isfinite(y)) {
			    mantissa_ = std::frexp(x*std::exp(y - M_LN2*(static_cast<int>(y/M_LN2) + 1)), &exponent_);
			    mantissa_ ? exponent_ += static_cast<int>(y/M_LN2) + 1 : exponent_ = std::numeric_limits<int>::min();
			} else            
			    throw(std::runtime_error("Zahl: argument of constructor is not a number"));
		};
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
		double toDouble() const { return std::ldexp(mantissa_, exponent_);};
		Zahl abs() { mantissa_ = std::abs(mantissa_); return *this;};
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
	
	inline Zahl pow(double arg) { 
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
}


#endif
