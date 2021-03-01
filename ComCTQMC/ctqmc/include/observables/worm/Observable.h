#ifndef CTQMC_INCLUDE_OBSERVABLES_WORM_OBSERVABLE_H
#define CTQMC_INCLUDE_OBSERVABLES_WORM_OBSERVABLE_H

#include "../../config/worm/Index.h"
#include "Measurement.h"

#include "../../Data.h"
#include "../../State.h"


namespace obs {
    
    namespace worm {
        

        template<typename Mode, typename Value, FuncType funcType, MeasType measType, typename Worm>
        struct Observable : obs::itf::Observable<Value> {
            Observable() = delete;
            Observable(std::int64_t store, jsx::value const& jWorm, data::Data<Value> const& data) :
            store_(store),
            jWorm_(jWorm),
            index_(data.ops().flavors()),
            samples0_(0),
            samples_(index_.size(),0),
            ptrs_(index_.size(), nullptr) {                                               
            };
            Observable(Observable const&) = delete;
            Observable(Observable&&) = delete;
            Observable& operator=(Observable const&) = delete;
            Observable& operator=(Observable&&) = delete;
            ~Observable() = default;
            
            bool sample(Value const sign, data::Data<Value> const& data, state::State<Value>& state, jsx::value& measurements, imp::itf::Batcher<Value>& batcher) {
                auto const integer_index = index_.integer(cfg::get<Worm>(state.worm()));
                
                if(ptrs_[integer_index] == nullptr) {
                    auto const string_index = cfg::worm::Index<Worm>::string(cfg::get<Worm>(state.worm()));
                    meas_.emplace(string_index, Meas<Mode, Value, funcType, measType, Worm>(jWorm_, data));
                    ptrs_[integer_index] = &meas_.at(string_index);
                    samplesMap_.emplace(string_index,&samples_[integer_index]);
                }
                
                ptrs_[integer_index]->add(sign, data, state);
                
                ++samples_[integer_index]; if(samples_[integer_index]%store_ == 0) store(data, measurements);
                
                return true;
            };
            
            void finalize(data::Data<Value> const& data, jsx::value& measurements) {
                std::size_t size = meas_.size(); std::size_t length = 0;
                for(auto const& basis : meas_) length = std::max(length, basis.first.size() + 1);
                mpi::all_reduce<mpi::op::max>(size); mpi::all_reduce<mpi::op::max>(length);
                
                std::vector<char> entries(size*length, '\0'); std::size_t pos = 0;
                for(auto const& basis : meas_) std::copy(basis.first.begin(), basis.first.end(), &entries[length*pos++]);
                mpi::gather(entries, mpi::master);
                
                std::vector<char> allEntries; std::int64_t allSize;
                if(mpi::rank() == mpi::master) {
                    std::set<std::string> guard;
                    for(pos = 0; pos*length < entries.size(); ++pos)
                        if(entries[pos*length] != '\0') guard.insert(std::string(&entries[pos*length]));
                    allEntries.resize(guard.size()*length, '\0'); pos = 0;
                    for(auto const& entry : guard) std::copy(entry.begin(), entry.end(), &allEntries[length*pos++]);
                    allSize = allEntries.size(); mpi::bcast(allSize, mpi::master);
                } else {
                    mpi::bcast(allSize, mpi::master); allEntries.resize(allSize);
                }
                mpi::bcast(allEntries, mpi::master);

                for(pos = 0; pos*length < allEntries.size(); ++pos) {
                    std::string const entry(&allEntries[pos*length]);
                    if(!meas_.count(entry)) meas_.emplace(entry, Meas<Mode, Value, funcType, measType, Worm>(jWorm_, data));
                }
                
                store(data, measurements);
            };
            
        private:
            std::int64_t const store_;
            jsx::value const jWorm_;
            
            cfg::worm::Index<Worm> index_;

            std::size_t samples0_;
            std::vector<std::int64_t> samples_;
            std::map<std::string,std::int64_t*> samplesMap_;
            
            std::vector<Meas<Mode, Value, funcType, measType, Worm>*> ptrs_;
            std::map<std::string, Meas<Mode, Value, funcType, measType, Worm>> meas_;
            
            
            void store(data::Data<Value> const& data, jsx::value& measurements) {
                
                for(auto& basis : meas_) {
                    auto& entry = measurements[measType == MeasType::Static ? "static" : "dynamic"][basis.first];
                    basis.second.store(entry, *samplesMap_.at(basis.first));
                    
                    *samplesMap_.at(basis.first)=0;
                }
                
            };
        };
        
    }
    
}

#endif
