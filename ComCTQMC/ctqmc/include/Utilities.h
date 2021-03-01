#ifndef CTQMC_INCLUDE_UTILITIES_H
#define CTQMC_INCLUDE_UTILITIES_H

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <complex>
#include <cstdint>
#include <random>
#include <limits>
#include <utility>
#include <stdexcept>


namespace ut {
    

    typedef std::complex<double> complex;
    
    inline double real(double arg) { return arg;};
    inline double real(complex arg) { return arg.real();};
    
    inline double conj(double arg) { return arg;};
    inline complex conj(complex arg) { return std::conj(arg);};
    
    template<typename T> inline T to_val(complex);
    
    template<> inline double to_val<double>(complex arg) { return arg.real();};
    template<> inline complex to_val<complex>(complex arg) { return arg;};
    
    
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
    
    enum class Flag { Pending, Accept, Reject };
    
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
    
    template<typename Value> struct Zahl;
    
    template<>
    struct Zahl<double> {
        Zahl() : mantissa_(.0), exponent_(std::numeric_limits<int>::min()) {};
        Zahl(double x, double y = .0) {   // constructs x*exp(y)
            if(std::isfinite(x) && std::isfinite(y)) {
                mantissa_ = std::frexp(x*std::exp(y - M_LN2*(static_cast<int>(y/M_LN2) + 1)), &exponent_);
                mantissa_ != .0 ? exponent_ += static_cast<int>(y/M_LN2) + 1 : exponent_ = std::numeric_limits<int>::min();
            } else
                throw std::runtime_error("ut::Zahl<double>: argument of constructor is not a number");
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
                mantissa_ != .0 ? exponent_ += exp : exponent_ = std::numeric_limits<int>::min();
            } else {
                mantissa_ = std::frexp(arg.mantissa_ + std::ldexp(mantissa_, exponent_ - arg.exponent_), &exp);
                mantissa_ != .0 ? exponent_ = arg.exponent_ + exp : exponent_ = std::numeric_limits<int>::min();
            }
            return *this;
        };
        Zahl& operator*=(Zahl const& arg) {
            int exp; mantissa_ = std::frexp(mantissa_*arg.mantissa_, &exp);
            mantissa_ != .0 ? exponent_ += (arg.exponent_ + exp) : exponent_ = std::numeric_limits<int>::min();
            return *this;
        };
        Zahl& operator/=(Zahl const& arg) {
            if(arg.mantissa_ == .0) throw std::runtime_error("ut::Zahl<double>: division by zero");
            
            int exp; mantissa_ = std::frexp(mantissa_/arg.mantissa_, &exp);
            mantissa_ != .0 ? exponent_ += (-arg.exponent_ + exp) : exponent_ = std::numeric_limits<int>::min();
            return *this;
        };
        double get() const {
            return std::ldexp(mantissa_, exponent_);
        };
        Zahl abs() const {
            Zahl temp = *this; temp.mantissa_ = std::abs(temp.mantissa_); return temp;
        };
        
        double mantissa() const { return mantissa_;};
        int exponent() const { return exponent_;};
    private:
        double mantissa_;
        int exponent_;
        
        friend struct Zahl<complex>;
        friend int operator==(Zahl<double> const&, Zahl<double> const&);
        friend int operator<=(Zahl<double> const&, Zahl<double> const&);
        friend Zahl<double> abs(Zahl<double> const&);
    };
    
    inline int operator==(Zahl<double> const& x, Zahl<double> const& y) {   // should this be implemented for complex case as well ?
        return (x.mantissa_ == y.mantissa_)&&(x.exponent_ == y.exponent_);
    }
    
    inline int operator<=(Zahl<double> const& x, Zahl<double> const& y) {
        if(x.mantissa_*y.mantissa_ <= .0 || x.exponent_ == y.exponent_) return x.mantissa_ <= y.mantissa_;
        return x.exponent_ < y.exponent_ ? x.mantissa_ > .0 : x.mantissa_ < .0;
    }
    
    inline Zahl<double> exp(double arg) {
        return arg != -std::numeric_limits<double>::infinity() ? Zahl<double>(1., arg) : Zahl<double>();
    }
    
    inline Zahl<double> operator+(Zahl<double> const& x, Zahl<double> const& y) {
        Zahl<double> temp(x); temp += y; return temp;
    }
    
    inline Zahl<double> operator*(Zahl<double> const& x, Zahl<double> const& y) {
        Zahl<double> temp(x); temp *= y; return temp;
    }
    
    inline Zahl<double> operator/(Zahl<double> const& x, Zahl<double> const& y) {
        Zahl<double> temp(x); temp /= y; return temp;
    }
    
    
    template<>
    struct Zahl<complex> {
        Zahl() : mantissa_(.0), exponent_(std::numeric_limits<int>::min()) {};
        Zahl(complex z, double y = .0) {
            if(std::isfinite(z.real()) && std::isfinite(z.imag()) && std::isfinite(y)) {  //better this way than using constructor of Zahl<double> because of throw
                double const x = std::abs(z);
                mantissa_ = x != .0 ? std::frexp(x*std::exp(y - M_LN2*(static_cast<int>(y/M_LN2) + 1)), &exponent_)*(z/x) : .0;
                mantissa_ != .0 ? exponent_ += static_cast<int>(y/M_LN2) + 1 : exponent_ = std::numeric_limits<int>::min();
            } else
                throw std::runtime_error("ut::Zahl<complex>: argument of constructor is not a number");
        };
        Zahl(Zahl const&) = default;
        Zahl(Zahl&&) = default;
        Zahl& operator=(Zahl const&) = default;
        Zahl& operator=(Zahl&&) = default;
        ~Zahl() = default;
        
        Zahl& operator+=(Zahl const& arg) {
            int exp;
            if(exponent_ > arg.exponent_) {
                auto const z = mantissa_ + (arg.mantissa_ != .0 ? arg.mantissa_*std::ldexp(1., arg.exponent_ - exponent_) : .0);
                auto const x = std::abs(z); mantissa_ = x != .0 ? std::frexp(x, &exp)*(z/x) : .0;
                mantissa_ != .0 ? exponent_ += exp : exponent_ = std::numeric_limits<int>::min();
            } else {
                auto const z = arg.mantissa_ + (mantissa_ != .0 ? mantissa_*std::ldexp(1., exponent_ - arg.exponent_) : .0);
                auto const x = std::abs(z); mantissa_ = x != .0 ? std::frexp(x, &exp)*(z/x) : .0;
                mantissa_ != .0 ? exponent_ = arg.exponent_ + exp : exponent_ = std::numeric_limits<int>::min();
            }
            return *this;
        };
        Zahl& operator*=(Zahl const& arg) {
            auto const z = mantissa_*arg.mantissa_; int exp;
            auto const x = std::abs(z); mantissa_ = x != .0 ? std::frexp(x, &exp)*(z/x) : .0;
            mantissa_ != .0 ? exponent_ += (arg.exponent_ + exp) : exponent_ = std::numeric_limits<int>::min();
            return *this;
        };
        Zahl& operator/=(Zahl const& arg) {
            if(arg.mantissa_ == .0) throw std::runtime_error("ut::Zahl<complex>: division by zero");
            
            auto const z = mantissa_/arg.mantissa_; int exp;
            auto const x = std::abs(z); mantissa_ = x != .0 ? std::frexp(x, &exp)*(z/x) : .0;
            mantissa_ != .0 ? exponent_ += (-arg.exponent_ + exp) : exponent_ = std::numeric_limits<int>::min();
            return *this;
        };
        complex get() const {
            return mantissa_*std::ldexp(1., exponent_);
        };
        Zahl<double> abs() const {
            Zahl<double> temp; temp.mantissa_ = std::abs(mantissa_); temp.exponent_ = exponent_; return temp;
        };
        
        complex mantissa() const { return mantissa_;};
        int exponent() const { return exponent_;};
    private:
        complex mantissa_;
        int exponent_;
        
        friend Zahl<double> abs(Zahl const&);
    };
    
    inline Zahl<complex> operator+(Zahl<complex> const& x, Zahl<complex> const& y) {
        Zahl<complex> temp(x); temp += y; return temp;
    }
    
    inline Zahl<complex> operator*(Zahl<complex> const& x, Zahl<complex> const& y) {
        Zahl<complex> temp(x); temp *= y; return temp;
    }
    
    inline Zahl<complex> operator/(Zahl<complex> const& x, Zahl<complex> const& y) {
        Zahl<complex> temp(x); temp /= y; return temp;
    }
    
    
    
    template<typename Value>
    inline Zahl<double> abs(Zahl<Value> const& arg) {
        return arg.abs();
    }
    
    
    //-------------------------------------------------------------------------------------------------------------------------
    
    template<typename...> using void_t = void;   // Walter Brown
    
    
    template<std::size_t...> struct sequence {};
    
    template<std::size_t Start, std::size_t End, std::size_t... Indices>   // sfinae for End < Start
    struct make_sequence {
        using type = typename make_sequence<Start, End - 1, End - 1, Indices...>::type;
    };
    
    template<std::size_t Start, std::size_t... Indices>
    struct make_sequence<Start, Start, Start, Indices...> {
        using type = sequence<Start, Indices...>;
    };
    
    template<std::size_t Start>
    struct make_sequence<Start, Start> : sequence<> {
        using type = sequence<>;
    };
    
    template<std::size_t Start, std::size_t End> using make_sequence_t = typename make_sequence<Start, End>::type;
    

    template<std::size_t, typename...> struct sequence_get;
    
    template<std::size_t Index, std::size_t Current, std::size_t... Rest>
    struct sequence_get<Index, sequence<Current, Rest...>> {
        constexpr static std::size_t value = sequence_get<Index - 1, sequence<Rest...>>::value;
    };
    
    template<std::size_t Current, std::size_t... Rest>
    struct sequence_get<0, sequence<Current, Rest...>> {
        constexpr static std::size_t value = Current;
    };
    
    
    template<typename...> struct transform_sequence;
    
    template<std::size_t... Input, std::size_t... Perm>
    struct transform_sequence<sequence<Input...>, sequence<Perm...>> {
        using type = sequence<sequence_get<Perm, sequence<Input...>>::value...>;
    };
    
    template<std::size_t... Input>
    struct transform_sequence<sequence<Input...>, sequence<>> {
        using type = sequence<>;
    };
    
    template<typename Input, typename Perm> using transform_sequence_t = typename transform_sequence<Input, Perm>::type;
    
    
    struct dummy_t {};  inline void sink(...) {}; // for using parameter pack expansion to iterate with std::get over std::tuple. Not very clean but handy ...
    
    
    template<typename T> struct wrap {};
    
    
    template<typename T, typename U, typename... Types>
    struct get_index_by_type {
        static constexpr int value = get_index_by_type<T, Types...>::value + 1;
    };
    
    template<typename T, typename... Types>
    struct get_index_by_type<T, T, Types...> {
        static constexpr int value = 0;
    };
    
    
    template<std::size_t Index, typename T, typename... Types>
    struct get_type_by_index {
        using type = typename get_type_by_index<Index - 1, Types...>::type;
    };
    
    template<typename T, typename... Types>
    struct get_type_by_index<0, T, Types...> {
        using type = T;
    };
    
    template<std::size_t Index, typename... Types> using get_type_by_index_t = typename get_type_by_index<Index, Types...>::type;
    
    
    template<typename...> struct for_each_type;
    
    template<typename T, typename... Types>
    struct for_each_type<T, Types...> {
        template<typename F, typename... Args>
        static void apply(F&& func, Args&&... args) {
            func(ut::wrap<T>(), std::forward<Args>(args)...);
            for_each_type<Types...>::apply(std::forward<F>(func), std::forward<Args>(args)...);
        }
    };
    
    template<>
    struct for_each_type<> {
        template<typename F, typename... Args>
        static void apply(F&& func, Args&&... fargs) {
        }
    };
    
    /*
    namespace impl {
        
        
        template<std::size_t Index, std::size_t Current, std::size_t Arg>
        struct replace_if {
            constexpr static std::size_t value = Arg;
        };
        
        template<std::size_t Index, std::size_t Current>
        struct replace_if<Index, Current, Index> {
            constexpr static std::size_t value = Current;
        };
        
        
        template<std::size_t...> struct check_if;
        
        template<std::size_t Index, std::size_t Current, std::size_t... Rest>
        struct check_if<Index, Current, Rest...> : check_if<Index, Rest...> {
        };
        
        template<std::size_t Index, std::size_t... Rest>
        struct check_if<Index, Index, Rest...> : std::true_type {
        };
        
        template<std::size_t Index>
        struct check_if<Index> : std::false_type {
        };
        
        
        template<std::size_t...> struct permutation_sign;
        
        template<std::size_t Index, std::size_t Current, std::size_t... Rest>
        struct permutation_sign<Index, Current, Rest...> {
            static_assert(check_if<Index, Rest...>::value, "not a permutation");
            
            static constexpr int value = -permutation_sign<Index + 1, replace_if<Index, Current, Rest>::value...>::value;
        };
        
        template<std::size_t Index, std::size_t... Rest>
        struct permutation_sign<Index, Index, Rest...> {
            static constexpr int value = permutation_sign<Index + 1, Rest...>::value;
        };
        
        template<std::size_t Index>
        struct permutation_sign<Index> {
            static constexpr int value = 1;
        };
        
    }
    
    
    template<std::size_t... Perm>
    struct permutation_sign {
        static constexpr int value = impl::permutation_sign<0, Perm...>::value;
    };
    */
}


#endif
