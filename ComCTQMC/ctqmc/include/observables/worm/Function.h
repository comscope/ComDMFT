#ifndef CTQMC_INCLUDE_OBSERVABLES_WORM_FUNCTION_H
#define CTQMC_INCLUDE_OBSERVABLES_WORM_FUNCTION_H


#include "../../Utilities.h"


namespace obs {
    
    namespace worm {

        
        enum class FuncType { Matsubara, Legendre };
        
        enum class PartType { Fermion, Boson };
        
        
        
        
        template<FuncType, PartType, typename> struct Function;
        

        template<PartType partType, typename Value>
        struct Function<FuncType::Legendre, partType, Value> {
            Function() = delete;
            Function(int size) : data_(size) {}
            Function(Function const&) = delete;
            Function(Function&&) = default;
            Function& operator=(Function const&) = delete;
            Function& operator=(Function&&) = delete;
            ~Function() = default;
    
            void operator()(ut::KeyType key) {
                
                bool antisymmetrize = false;
                if(key < 0) {
                    key += ut::KeyMax;
                    antisymmetrize = true;
                }
                
                data_[0] = 1.;
                data_[1] = ((2.*key)/ut::KeyMax - 1.);
                
                double const x = data_[1];
                for(std::size_t l = 2; l < data_.size(); ++l)
                    data_[l] = ((2*l - 1)*x*data_[l - 1] - (l - 1)*data_[l - 2])/l;
                
                if (antisymmetrize){
                    constexpr double p = partType == PartType::Fermion ? -1. : 1.;
                    for (auto& x : data_) x*=p;
                }
            };
            
            std::vector<double> const& operator()() {
                return data_;
            };
            
        private:
            std::vector<double> data_;
        };
        
        
        
        
        template<typename Iterator>
        void set_matsubara(Iterator start, Iterator end, ut::complex val, ut::complex const fact) {
            for(Iterator it = start; it != end; ++it) {
                *it = val; val *= fact;
            }
        }
        
        
        template<PartType partType>
        struct Function<FuncType::Matsubara, partType, double> {
            Function() = delete;
            Function(int size) : data_(size) {};
            Function(Function const&) = delete;
            Function(Function&&) = default;
            Function& operator=(Function const&) = delete;
            Function& operator=(Function&&) = delete;
            ~Function() = default;
            
            void operator()(ut::KeyType key) {
                double const phi = (M_PI*key)/ut::KeyMax;
                
                ut::complex const fact{ std::cos(2.*phi), std::sin(2.*phi) };
                ut::complex const val{ partType == PartType::Fermion ? ut::complex{ std::cos(phi), std::sin(phi) } : 1. };

                set_matsubara(data_.begin(), data_.end(), val, fact);
            };
            
            std::vector<ut::complex> const& operator()() {
                return data_;
            };
            
        private:
            std::vector<ut::complex> data_;
        };
        

        
        
        template<PartType partType>
        struct Function<FuncType::Matsubara, partType, ut::complex> {
            Function() = delete;
            Function(int size) : data_(partType == PartType::Fermion ? 2*size : 2*size - 1) {}; // size is the number of positive matsubara frequencies (including 0 if bosonic)
            Function(Function const&) = delete;
            Function(Function&&) = default;
            Function& operator=(Function const&) = delete;
            Function& operator=(Function&&) = delete;
            ~Function() = default;
            
            void operator()(ut::KeyType key) {
                double const phi = (M_PI*key)/ut::KeyMax;
                
                ut::complex const fact{ std::cos(2.*phi), std::sin(2.*phi) };
                ut::complex const val{ partType == PartType::Fermion ? ut::complex{ std::cos(phi), std::sin(phi) } : 1. };
                std::size_t const start = data_.size()/2; // start = size if fermion and size - 1 if boson, so we are setting bosonic frequency zero twice ... who cares.

                set_matsubara( data_.begin() + start,  data_.end(),           val ,           fact );
                set_matsubara(data_.rbegin() + start, data_.rend(), std::conj(val), std::conj(fact));  
            };
            
            std::vector<ut::complex> const& operator()() {
                return data_;
            };
            
        private:
            std::vector<ut::complex> data_;
        };

        
        
    }
    
}

#endif
