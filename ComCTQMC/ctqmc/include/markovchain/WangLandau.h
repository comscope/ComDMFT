#ifndef CTQMC_INCLUDE_MARKOVCHAIN_WANGLANDAU_H
#define CTQMC_INCLUDE_MARKOVCHAIN_WANGLANDAU_H


#include "../Data.h"
#include "../config/Worms.h"
#include "../config/worm/Index.h"
#include "../../../include/JsonX.h"

namespace mch {
    
    
    template<typename Value>
    struct WangLandau {
        WangLandau() = delete;
        WangLandau(jsx::value const& jParams, data::Data<Value> const& data) :
        restart_(jParams.is("restart") and jParams("restart").boolean()),
        names_(cfg::Worm::get_names()),
        flatCriterion_(0.4),
        zFrac_(jParams.is("partition fraction") ?
                                             jParams("partition fraction").real64() : 0.5),
        number_of_spaces_(0),
        lambda_(2.),
        thermalised_(false), totalSteps_(0) {
            steps_.resize(cfg::Worm::size()); eta_.resize(cfg::Worm::size());
            
            int index = 0;
            for(auto const& name : names_) {
                if(jParams.is(name)) {
                    active_.push_back(index);
                    
                    if (restart_){
                        auto const& jWangLandau = jParams("measurements")("WL");
                        if(jWangLandau.is(name)){
                        
                            for (auto const& worm : jWangLandau(name).object())
                                safe_emplace(index, worm.first,
                                                restart_ ? worm.second("eta").real64() : 1.,
                                                restart_ ? worm.second("steps").int64()/mpi::number_of_workers() : 0);
                            
                        }
                    }

                }
                ++index;
            }
             
            auto const it = std::find(active_.begin(), active_.end(), index = get_index<cfg::partition::Worm>::value);
            if( it == active_.end())
                throw std::runtime_error("mch::WangLandau: partition space not active !");
             
        };
        WangLandau(WangLandau const&) = delete;
        WangLandau(WangLandau&&) = delete;
        WangLandau& operator=(WangLandau const&) = delete;
        WangLandau& operator=(WangLandau&&) = delete;
        ~WangLandau() = default;
        
        void thermalised() {
            
            //Guard against calling this multiple times -- which happens when there are GPU
            if (!thermalised_){
                    
                if(!restart_){
                    normalise_eta();  fill<std::int64_t>(steps_,0); totalSteps_ = 0;
                }
                
                //All mp images should have the same eta
                for(auto active : active_)
                    for (auto & eta : eta_[active]){
                        mpi::reduce<mpi::op::sum>(eta.second, mpi::master);
                    
                        if(mpi::rank() == mpi::master)
                            eta.second /= mpi::number_of_workers();
                        
                        mpi::bcast(eta.second, mpi::master);
                    }
                
                
                // let user decide how much time should be spent in partition space (default 50%)
                /* eqn follows frmo starts from
                 \eta_z' V_z / (\eta_z' V_z + \Sum_O \eta_O V_O) = zFrac_       (1)
                 and
                 \eta_z V_z / (\eta_z V_z + \Sum_O \eta_O V_O = 1 / N_{spaces}  (2)
                 
                 // from (2) input \eta V_z into (1)
                */
                
                if (!restart_ and number_of_spaces_ > 1)
                    eta_[get_index<cfg::partition::Worm>::value].at(cfg::worm::Index<cfg::partition::Worm>::string()) *= (number_of_spaces_ - 1.) / (1./zFrac_ - 1.);
                
                double const partition = eta_[get_index<cfg::partition::Worm>::value].at(cfg::worm::Index<cfg::partition::Worm>::string());

                //Let user know what eta's were chosen
                //and normalize to eta_z=1
                for(auto active : active_){
                    
                    mpi::cout << names_[active] << ": ";
                    
                    for (auto & eta : eta_[active]){
                        eta.second /= partition;  // Tradition ...
                    
                        mpi::cout << "eta_" << eta.first << " = " << eta.second << " | ";
                    }
                    
                    mpi::cout << std::endl;
                }
                
                thermalised_ = true;
            }
        };
        
        template<typename W>
        double eta(W const& w) const {
            auto const it = eta_[get_index<W>::value].find(cfg::worm::Index<W>::string(w));
            return it != eta_[get_index<W>::value].end() ? it->second : 1.;
        };
        
        template<typename W>
        std::int64_t steps(std::string const& string_index) const {
            auto const it = steps_[get_index<W>::value].find(string_index);
            return it != steps_[get_index<W>::value].end() ? it->second : 0;
        };
        
        template<typename W>
        void inc(W const& w) {
            
            auto const string_index = cfg::worm::Index<W>::string(w);
            safe_emplace(get_index<W>::value,string_index);
            
            if (!restart_ or thermalised_){ ++steps_[get_index<W>::value].at(string_index); ++totalSteps_; }
            
            if (!thermalised_ and !restart_) {
                
                //Update eta, avoid underflows ...
                if(eta_[get_index<W>::value].at(string_index) > 1.e-25) {
                    eta_[get_index<W>::value].at(string_index) /= lambda_;
                    //for (auto& eta : eta_[get_index<W>::value])
                    //    eta.second = eta_[get_index<W>::value].at(string_index);
                    
                } else {
                    normalise_eta();
                    
                    if(eta_[get_index<W>::value].at(string_index) > 1.e-25) {
                        eta_[get_index<W>::value].at(string_index) /= lambda_;
                        //for (auto& eta : eta_[get_index<W>::value])
                            //eta.second = eta_[get_index<W>::value].at(string_index);
                        
                    } else {
                        lambda_ = std::sqrt(lambda_);  fill(eta_,1.); fill<std::int64_t>(steps_,0); totalSteps_=0;  return; // ... and restart with smaller lambda if eta's get crazy
                    }
                }
                
                //Look at the current histogram ...
                std::int64_t min = std::numeric_limits<std::int64_t>::max(), max = 0;
                for(auto const active : active_) {
                    if (steps_[active].size()){
                        for (auto & steps : steps_[active]){
                            min = std::min(steps.second, min);  max = std::max(steps.second, max);
                        }
                    } else {
                        min = 0; //Needed for startup
                    }
                }
                
                // ... is it sufficiently flat?  if so, update lambda and reset histogram
                if(max - min < totalSteps_*flatCriterion_/number_of_spaces_) {
                    lambda_ = std::sqrt(lambda_);  fill<std::int64_t>(steps_,0);  totalSteps_ = 0;
                }
            }
            
        };
        
        jsx::value json() {
            jsx::value jWangLandau;
            
            for(auto active : active_) {
                for (auto & steps : steps_[active]){
                    
                    jWangLandau[names_[active]][steps.first]["steps"] = steps.second;
                    
                    jWangLandau[names_[active]][steps.first]["eta"] = eta_[active].at(steps.first);
                    
                }
            }
            
            return jWangLandau;
        };
        
        
    private:
        template<typename W> using get_index = cfg::get_index<W, cfg::Worm>;
        bool const restart_;
        
        std::vector<std::string> const names_;
        double const flatCriterion_;  // Standard deviation at which point lambda is updated
        double const zFrac_;  // fraction of time one should spend in partition space
        std::size_t number_of_spaces_; //total number of configuratioon spaces
        
        double lambda_;               //Update the volume as  V -> V*lambda
        bool thermalised_;
        std::int64_t totalSteps_;
        
        std::vector<std::map<std::string,double>> eta_;
        std::vector<std::map<std::string,std::int64_t>> steps_;
        
        std::vector<int> active_;
        
        void safe_emplace(std::size_t const worm_index, std::string const& worm_entry, double const eta = 1., std::int64_t const steps = 0){
            {
                auto const it = eta_[worm_index].find(worm_entry);
                if (it == eta_[worm_index].end()){
                    eta_[worm_index].emplace(worm_entry,eta);
                }
            }
            {
                auto const it = steps_[worm_index].find(worm_entry);
                if (it == steps_[worm_index].end()){
                    steps_[worm_index].emplace(worm_entry,steps);
                    ++number_of_spaces_;
                }
            }
           
        }
        
        template <typename T>
        void fill(std::vector<std::map<std::string,T>> & vector_of_maps, T const arg ) {
            
            for (auto & map : vector_of_maps)
                for (auto & el : map)
                    el.second=arg;
            
        }
        
        void normalise_eta() {
            double norm = .0;
            for (auto const& eta_map : eta_)
                for (auto const& eta : eta_map)
                    norm += eta.second*eta.second;
            
            norm = std::sqrt(norm);
            for (auto & eta_map : eta_)
                for (auto & eta : eta_map){
                    eta.second /= norm;
                }
        }
    };

}





/*
 Index<Worm> index_;
 Index<Worm>::string(w)
 auto const integer_index = index_.integer(cfg::get<Worm>(state.worm()));
 */

#endif
