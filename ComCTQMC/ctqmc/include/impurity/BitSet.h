#ifndef CTQMC_INCLUDE_IMPURITY_BITSET_H
#define CTQMC_INCLUDE_IMPURITY_BITSET_H

#include <vector>

namespace imp {

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
        ~BitSet() { delete[] data_;};
        
        int any() const {
            for(std::size_t i = 0; i < size_; ++i) if(data_[i]) return 1;
            return 0;
        };
        int operator[](std::size_t pos) const {
            return data_[pos/sizeof(data_type)]&(static_cast<data_type>(1) << pos%sizeof(data_type));
        };
        
        void set(std::size_t pos) {
            data_[pos/sizeof(data_type)] |= (static_cast<data_type>(1) << pos%sizeof(data_type));
        };
       
    private:
        std::size_t const size_;
        data_type* const data_;
    };
}


#endif
