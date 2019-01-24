#ifndef UTILITIES
#define UTILITIES

#include <vector>
#include <complex>
#include <fstream>
#include <valarray>
#include <cstring>
#include <cstdint>
#include <random>

//Achtung: es kann sein dass gewisse observabeln nicht gespeichert wurden, c.f. MonteCarlo.h

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
        
        double operator()() const { return *beta_;};
        
        ~Beta() {
            delete beta_;
        };
    private:
        double* beta_;
    };
    
    Beta beta;
    
    //--------------------------------------------------------------------------------------------------------------------------

    struct BitSet {
        typedef std::uint_fast64_t data_type;
        
        BitSet() = delete;
        BitSet(std::size_t n) : size_((n + sizeof(data_type) - 1)/sizeof(data_type)), data_(new data_type[size_]) {
            std::memset(data_, 0, size_*sizeof(data_type));
        };
        BitSet(BitSet const&) = delete;
        BitSet(BitSet&&) = delete;
        BitSet& operator=(BitSet const&) = delete;
        BitSet& operator=(BitSet&&) = delete;
        
        int any() {
            for(std::size_t i = 0; i < size_; ++i) if(data_[i]) return 1;
            return 0;
        };
        int operator[](std::size_t pos) {
            return data_[pos/sizeof(data_type)]&(static_cast<data_type>(1) << pos%sizeof(data_type));
        };
        void flip(std::size_t pos) {
            data_[pos/sizeof(data_type)] ^= (static_cast<data_type>(1) << pos%sizeof(data_type));
        };
        ~BitSet() { delete[] data_;};
    private:
        std::size_t size_;
        data_type* data_;
    };
    
}


#endif
