#ifndef INCLUDE_IO_TENSOR_H
#define INCLUDE_IO_TENSOR_H

#include <vector>
#include <complex>

#include "Vector.h"
#include "../JsonX.h"
#include "../../ctqmc/include/Utilities.h"

//scheiss b64 member variable !! Kack loesig im moment ...

namespace io {
    
    template<typename T> struct Tensor {
        inline static std::string name() { return name(T());};
        
        Tensor() = default;
        Tensor(std::size_t I, std::size_t J, std::size_t K, std::size_t L) : I_(I), J_(J), K_(K), L_(L) {};
        Tensor(Tensor const&) = default;
        Tensor(Tensor&& other) noexcept : I_(other.I_), J_(other.J_),K_(other.K_), L_(other.L_), ijkl_(std::move(other.ijkl_)), data_(std::move(other.data_)), entries_(std::move(other.entries_)) { other.I_ = other.J_ = other.K_ = other.L_ = 0;};
        Tensor& operator=(Tensor const&) = default;
        Tensor& operator=(Tensor&& other) { I_ = other.I_; J_ = other.J_; K_ = other.K_; L_ = other.L_; other.I_ = other.J_ = other.K_ = other.L_ = 0; data_ = std::move(other.data_); entries_ = std::move(other.entries_); ijkl_ = std::move(other.ijkl_);  return *this;};
        ~Tensor() = default;

        int const& I() const { return I_;};
        int const& J() const { return J_;};
        int const& K() const { return K_;};
        int const& L() const { return L_;};
        
        T* data() { return data_.data();};
        T const* data() const { return data_.data();};
        std::vector<std::vector<int>> const& ijkl() const { return ijkl_;}
        
        
        T& operator()(int const i, int const j, int const k, int const l){ return this->operator()(this->index(i,j,k,l)); }
        T const& operator()(int const i, int const j, int const k, int const l) const { return this->operator()(this->index(i,j,k,l)); }
        
        T const at(int const i, int const j, int const k, int const l) const { return this->at(this->index(i,j,k,l)); }
        
        std::string const& entry(int const i, int const j, int const k, int const l) const { return this->entry(this->index(i,j,k,l)); }
        bool const& is(int const i, int const j, int const k, int const l) const { return this->is(this->index(i,j,k,l)); }
        
        
        inline T& operator()(int const i){
            if (data_.find(i) == data_.end())
                throw std::runtime_error("Cannot update element of " + this->name() + " that does not yet exist\n");
            return data_.at(i);
        };
        inline T const& operator()(int i) const {
            if (data_.find(i) != data_.end())
                return data_.at(i);
            else
                return zero;
        };
        
        inline T const at(int i) const {
            if (data_.find(i) != data_.end())
                return data_.at(i);
            else
                return zero;
        };
        
        inline std::string const& entry(int i) const {
            if (entries_.find(i) != entries_.end())
                return entries_.at(i);
            else
                return empty_entry;
        };
                              

        inline bool is(int i) const { return (data_.find(i) == data_.end()) ? false : true; }

        inline void emplace (int const i, int const j, int const k, int const l, std::string const& entry, T const v){
            
            ijkl_.push_back(std::vector<int>({this->index(i,j,k,l),i,j,k,l}));
            this->emplace(this->index(i,j,k,l),entry,v);
        }
        
        Tensor& resize(int I, int J, int K, int L, T value = .0) {
            return *this;
        };
        Tensor& conj() {
            Tensor temp;
            temp.I_ = J_; temp.J_ = I_; temp.K_ = K_; temp.L_ = L_;
            
            for(int i = 0; i < I_; ++i)
                for(int j = 0; j < J_; ++j)
                    for(int k = 0; k < K_; ++k)
                        for(int l = 0; l < L_; ++l){
                            temp(l, k, j, i) = ut::conj(operator()(i, j, k, l));
                        }
            
            return *this = std::move(temp);
        };
        
        /*
        void read(jsx::value const& source) {
            I_ = source(0).int64();
            J_ = source(1).int64();
            K_ = source(2).int64();
            L_ = source(3).int64();
            data_.read(source(4));
        };
        void write(jsx::value& dest) const {
            dest = jsx::array_t(5);
            
            dest(0) = static_cast<jsx::int64_t>(I());
            dest(1) = static_cast<jsx::int64_t>(J());
            dest(2) = static_cast<jsx::int64_t>(K());
            dest(3) = static_cast<jsx::int64_t>(L());
            data_.write(dest(4));
        };
        bool& b64() const {
            return data_.b64();
        };*/
        
    private:
        int I_ = 0, J_ = 0, K_ = 0, L_ = 0;
        std::vector<std::vector<int>> ijkl_;
        std::map<int,T> data_;
        std::map<int,std::string> entries_;

        inline int index(int const i, int const j, int const k, int const l) const { return i + j*I_ + k*I_*J_ + l*I_*J_*K_;}
        inline void emplace (int const i, std::string const& entry, T const v){
            entries_.emplace(i,entry);
            data_.emplace(i,v);
        }
        
        T const zero = 0.;
        std::string const empty_entry = "";
        
        inline static std::string name(double const&) { return "io::rtens";};
        inline static std::string name(std::complex<double> const&) { return "io::ctens";};
        
    };
    
    
    typedef Tensor<double> rtens;
    typedef Tensor<std::complex<double>> ctens;
    
    
};

#endif 
