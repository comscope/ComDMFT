#ifndef CTQMC_INCLUDE_OBSERVABLES_WORM_BASIS_H
#define CTQMC_INCLUDE_OBSERVABLES_WORM_BASIS_H

#include "Function.h"
#include "../../Data.h"


namespace obs {
    
    namespace worm {
        
        
        template<FuncType, typename Value> struct data_trait;
        
        template<typename Value>
        struct data_trait<FuncType::Matsubara, Value> {
            using type = ut::complex;
        };
        
        template<typename Value>
        struct data_trait<FuncType::Legendre, Value> {
            using type = Value;
        };
        
        template<FuncType funcType, typename Value> using data_trait_t = typename data_trait<funcType, Value>::type;
        
        
        
        
        template<typename, FuncType, typename...> struct Basis;
        
        template<typename Value, FuncType funcType>
        struct Basis<Value, funcType, cfg::FermionicTime, cfg::FermionicTime> {

            void store(jsx::value& measurements, std::int64_t samples) {
                measurements << meas::fix(data_, samples);
                
                std::fill(data_.begin(), data_.end(), .0);
            };
            
        protected:
            Basis() = delete;
            Basis(jsx::value const& jWorm, data::Data<Value> const& data) :
            fermion_(jWorm("cutoff").int64()),
            data_(fermion_().size(), .0) {
            };
            Basis(Basis const&) = delete;
            Basis(Basis&&) = default;
            Basis& operator=(Basis const&) = delete;
            Basis& operator=(Basis&&) = delete;
            ~Basis() = default;
            
            void add(Value const sign, cfg::FermionicTime const& op1, cfg::FermionicTime const& op2) {
                fermion_(op1.key() - op2.key());
                
                auto data = data_.data();
                for(auto f : fermion_())
                    *data++ += sign*f;
            };
            
            Function<funcType, PartType::Fermion, Value> fermion_;
            
            std::vector<data_trait_t<funcType, Value>> data_; // we need a type trait here to decide whether data_ is real or complex
        };
        
        
        template<typename Value, FuncType funcType>
        struct Basis<Value, funcType, cfg::FermionicTime, cfg::FermionicTime, cfg::FermionicTime, cfg::FermionicTime> {
            
            void store(jsx::value& measurements, std::int64_t samples) {
                measurements << meas::fix(data_, samples);
                
                std::fill(data_.begin(), data_.end(), .0);
            };
            
        protected:
            Basis() = delete;
            Basis(jsx::value const& jWorm, data::Data<Value> const& data) :
            boson32_(jWorm("boson cutoff").int64()),
            fermion12_(jWorm("fermion cutoff").int64()),
            fermion34_(jWorm("fermion cutoff").int64()),
            data_(fermion12_().size()*fermion34_().size()*boson32_().size(), .0) {
            };
            Basis(Basis const&) = delete;
            Basis(Basis&&) = default;
            Basis& operator=(Basis const&) = delete;
            Basis& operator=(Basis&&) = delete;
            ~Basis() = default;
            
            void add(Value const sign, cfg::FermionicTime const& op1, cfg::FermionicTime const& op2, cfg::FermionicTime const& op3, cfg::FermionicTime const& op4) {
                fermion12_(op1.key() - op2.key());
                fermion34_(op3.key() - op4.key());
                boson32_(op2.key() - op3.key());
                
                auto data = data_.data();
                for(auto b32 : boson32_())
                    for(auto f34 : fermion34_())
                        for(auto f12 : fermion12_())
                            *data++ += sign*f12*f34*b32;
            };
            
            Function<FuncType::Matsubara,   PartType::Boson,       Value> boson32_;
            Function<           funcType, PartType::Fermion, ut::complex> fermion12_;
            Function<           funcType, PartType::Fermion, ut::complex> fermion34_;

            std::vector<ut::complex> data_;   // always complex
        };
        
        
        template<typename Value, FuncType funcType>
        struct Basis<Value, funcType, cfg::BosonicTime, cfg::BosonicTime> {    // Fermionic and Bosonic one time class could be derived from parent class, but I do NOT think that this improves readability
            
            void store(jsx::value& measurements, std::int64_t samples) {
                measurements << meas::fix(data_, samples);
                
                std::fill(data_.begin(), data_.end(), .0);
            };
            
        protected:
            Basis() = delete;
            Basis(jsx::value const& jWorm, data::Data<Value> const& data) :
            boson_(jWorm("cutoff").int64()),
            data_(boson_().size(), .0) {
            };
            Basis(Basis const&) = delete;
            Basis(Basis&&) = default;
            Basis& operator=(Basis const&) = delete;
            Basis& operator=(Basis&&) = delete;
            ~Basis() = default;
            
            void add(Value const sign, cfg::BosonicTime const& op1, cfg::BosonicTime const& op2) {
                boson_(op1.key() - op2.key());
                
                auto data = data_.data();
                for(auto b : boson_())
                    *data++ += sign*b;
            };
            
            Function<funcType, PartType::Boson, Value> boson_;
            
            std::vector<data_trait_t<funcType, Value>> data_;  // we need a type trait here to decide whether data_ is real or complex
        };
        
        
        template<typename Value, FuncType funcType>
        struct Basis<Value, funcType, cfg::FermionicTime, cfg::FermionicTime, cfg::BosonicTime> {
            
            void store(jsx::value& measurements, std::int64_t samples) {
                measurements << meas::fix(data_, samples);
                
                std::fill(data_.begin(), data_.end(), .0);
            };
            
        protected:
            Basis() = delete;
            Basis(jsx::value const& jWorm, data::Data<Value> const& data) :
            boson_(jWorm("boson cutoff").int64()),
            fermion_(jWorm("fermion cutoff").int64()),
            data_(fermion_().size()*boson_().size(), .0){
            };
            Basis(Basis const&) = delete;
            Basis(Basis&&) = default;
            Basis& operator=(Basis const&) = delete;
            Basis& operator=(Basis&&) = delete;
            ~Basis() = default;
            
            void add(Value const sign, cfg::FermionicTime const& op1, cfg::FermionicTime const& op2, cfg::BosonicTime const& op3) {
                fermion_(op1.key() - op2.key());
                boson_(op2.key() - op3.key());
                
                auto data = data_.data();
                for(auto b : boson_())
                    for(auto f : fermion_())
                        *data++ += sign*f*b;
            };

            Function<FuncType::Matsubara,   PartType::Boson,       Value> boson_;
            Function<           funcType, PartType::Fermion, ut::complex> fermion_;  // I changed this back to ut::complex since we always need the positive and negative matsubara frequencies

            std::vector<ut::complex> data_;   // always complex
        };
        
    }
    
}

#endif
