#ifndef INCLUDE_MEASUREMENTS_MEASUREMENTS_H
#define INCLUDE_MEASUREMENTS_MEASUREMENTS_H

#include <vector>
#include <complex>
#include <fstream>
#include <valarray>
#include <cstring>
#include <random>

#include "../JsonX.h"
#include "../mpi/Utilities.h"
#include "../io/Vector.h"
#include "../../ctqmc/include/config/Worms.h"
#include "../../ctqmc/include/config/worm/Index.h"

//Achtung: es kann sein dass gewisse observabeln nicht gespeichert wurden, c.f. MonteCarlo.h

namespace meas {
    
    struct All {}; struct Jackknife {};
    
    
    struct Fix {};  struct Var {};

    template<typename T, typename M>
    struct Vector {
        inline static std::string name() { return name(T(), M());};
        
        Vector() = default;
        
        /*
        Vector(io::Vector<T> && data, std::int64_t const samples) : data_(0), samples_(0) {
            if(mpi::rank() == mpi::master){
                data_ = std::move(data);
                samples_ = samples;
            }
        }
        */
        
        Vector(Vector const&) = default;
        Vector(Vector&&) = default;
        Vector& operator=(Vector const&) = default;
        Vector& operator=(Vector&&) = default;
        ~Vector() = default;
        
        void add(std::vector<T> const& val, std::int64_t samples) {
            resize_add(val, data_, M()); samples_ += samples; for(std::size_t n = 0; n < val.size(); ++n) data_[n] += val[n];
        };
        
        jsx::value reduce(double fact, All, bool b64) const {
            auto samples = samples_;  mpi::reduce<mpi::op::sum>(samples, mpi::master);
            auto data = data_;  resize_reduce(data, M());  mpi::reduce<mpi::op::sum>(data, mpi::master);
            
            if(mpi::rank() == mpi::master) {
                if(!samples) throw std::runtime_error(name() + "::write: no measurements taken !");  //Scheisse das sštt eh nie passiere.
                for(auto& x : data) x *= fact/samples;
                data.b64() = b64;  return std::move(data);
            } else
                return jsx::null_t();
        };
        
        jsx::value reduce(double fact, Jackknife, bool b64) const {
            auto samples = samples_;  mpi::all_reduce<mpi::op::sum>(samples);  samples -= samples_;
            
            if(!samples) throw std::runtime_error(name() + "::write: no measurements taken !");  //Scheisse das sštt eh nie passiere.
            
            auto data = data_;  resize_reduce(data, M());  mpi::all_reduce<mpi::op::sum>(data);
            for(std::size_t i = 0; i < data_.size(); ++i) data[i] -= data_[i];
            for(auto& x : data) x *= fact/samples;
            
            data.b64() = b64; return std::move(data);
        };
        
        void write(jsx::value& dest) const {
            throw std::runtime_error(name() + "::write: not implemented");
        };
        
    private:
        
        std::int64_t samples_ = 0;
        io::Vector<T> data_;
        
        static std::string name(double const&, Fix)               { return "meas::rvecfix"; };
        static std::string name(std::complex<double> const&, Fix) { return "meas::cvecfix"; };
        
        static std::string name(double const&, Var)               { return "meas::rvecvar"; };
        static std::string name(std::complex<double> const&, Var) { return "meas::cvecvar"; };
        
        static void resize_add(io::Vector<T> const& val, io::Vector<T>& data, Fix) {
            if(data.size() == 0) data.resize(val.size(), .0);
            if(data.size() != val.size()) throw std::runtime_error(name() + "::add: missmatch in array size!");
        };
        static void resize_add(io::Vector<T> const& val, io::Vector<T>& data, Var) {
            if(val.size() > data.size()) data.resize(val.size(), .0);
        };
        
        static void resize_reduce(io::Vector<T>& data, Fix) {
        };
        static void resize_reduce(io::Vector<T>& data, Var) {
            auto size = data.size(); mpi::all_reduce<mpi::op::max>(size); data.resize(size, .0);
        };
    };

    
    template<typename T, typename M>
    struct Sample {
        T const& value; std::int64_t samples;
    };
    
    
    template<typename T> inline Sample<T, Fix> fix(T const& value, std::int64_t samples) { return {value, samples};};
    template<typename T> inline Sample<T, Var> var(T const& value, std::int64_t samples) { return {value, samples};};
    
    
    template<typename T, typename M>
    inline void operator<<(jsx::value& lhs, Sample<T, M> const& rhs) {
        if(lhs.is<jsx::empty_t>()) lhs = Vector<T, M>();
        lhs.at<Vector<T, M>>().add(std::vector<T>(1, rhs.value), rhs.samples);
    }
    
    template<typename T, typename M>
    inline void operator<<(jsx::value& lhs, Sample<std::vector<T>, M> const& rhs) {
        if(lhs.is<jsx::empty_t>()) lhs = Vector<T, M>();
        lhs.at<Vector<T, M>>().add(rhs.value, rhs.samples);
    }

    
    using rvecfix = Vector<double, Fix>;  using cvecfix = Vector<std::complex<double>, Fix>;
    using rvecvar = Vector<double, Var>;  using cvecvar = Vector<std::complex<double>, Var>;
    
    
    //--------------------------------------------------------------------------------------------------------------------------------
    
    
    std::int64_t reduce_steps(std::int64_t steps, All) {
        mpi::reduce<mpi::op::sum>(steps, mpi::master); return steps;
    };
    
    std::int64_t reduce_steps(std::int64_t steps, Jackknife) {
        auto temp = steps; mpi::all_reduce<mpi::op::sum>(temp); return temp - steps;
    };
    
    
    template<typename E>
    inline void reduce(jsx::value& jOut, double fact, jsx::value const& jIn, E, bool b64) {
        if(jIn.is<rvecfix>())
            jOut = jIn.at<rvecfix>().reduce(fact, E(), b64);
        else if(jIn.is<cvecfix>())
            jOut = jIn.at<cvecfix>().reduce(fact, E(), b64);
        else if(jIn.is<rvecvar>())
            jOut = jIn.at<rvecvar>().reduce(fact, E(), b64);
        else if(jIn.is<cvecvar>())
            jOut = jIn.at<cvecvar>().reduce(fact, E(), b64);
        else if(jIn.is<jsx::object_t>()) {
            for(auto& jEntry : jIn.object()) reduce(jOut[jEntry.first], fact, jEntry.second, E(), b64);
        } else if(jIn.is<jsx::array_t>()) {
            if(!(jOut.is<jsx::array_t>() && jOut.size() == jIn.size())) jOut = jsx::array_t(jIn.size());
            int index = 0; for(auto& jEntry : jIn.array()) reduce(jOut[index++], fact, jEntry, E(), b64);
        } else
            jOut = jIn;
    }
    
    
    template<typename E>
    inline void reduce(jsx::value& jOut, jsx::value const& jIn, jsx::value const& jWL, E, bool b64) {
        auto const pName = cfg::partition::Worm::name();
        auto const pSpec = cfg::worm::Index<cfg::partition::Worm>::string();
        
        jsx::value const jSign = jIn(pName)("sign").at<rvecfix>().reduce(1., E(), b64);
        auto const pSteps = reduce_steps(jWL(pName)(pSpec)("steps").int64(), E());
        auto const pEta = jWL(pName)(pSpec)("eta").real64();
        auto const signxZp = (jSign.is<jsx::null_t>() ? 1. : jsx::at<io::rvec>(jSign).at(0))*pSteps/pEta;
        
        if(!jOut.is<jsx::object_t>()) jOut = jsx::object_t();
        
        for(auto& jWorm : jIn.object()) {
            if (jWorm.first==pName){
                reduce(jOut[jWorm.first], pSteps/pEta/signxZp, jWorm.second, E(), b64);
                jOut["WL"][jWorm.first][pSpec]["steps"] = pSteps;
                jOut["WL"][jWorm.first][pSpec]["eta"] = pEta;
            } else {
                for (auto& type : jWorm.second.object())
                    for (auto& func : type.second.object()){
                        auto const wSteps = reduce_steps(jWL(jWorm.first)(func.first)("steps").int64(), E());
                        auto const wEta = jWL(jWorm.first)(func.first)("eta").real64();
                        auto const Zw = wSteps/wEta;
                            
                        reduce(jOut[jWorm.first][type.first][func.first], Zw/signxZp, func.second, E(), b64);
                        jOut["WL"][jWorm.first][func.first]["steps"] = wSteps;
                        jOut["WL"][jWorm.first][func.first]["eta"] = wEta;
                    }
            }
        }
        
        jOut[pName]["sign"] = jSign;
    }

    void restart(jsx::value const& jParams, jsx::value & jInput, jsx::value & jMeasurements){
        if(mpi::rank() == mpi::master){
            mpi::cout << "Initializing measurements from previous run" << std::endl;
            auto const names = cfg::Worm::get_names();
            auto const pName = cfg::partition::Worm::name();
            auto const pSpec = cfg::worm::Index<cfg::partition::Worm>::string();
            
            for (auto & meas : jInput.object()){
                
                //Only handle the actual measurements in this way.
                //Other information handled separately, e.g., WangLandau steps and etas
                auto it = std::find(names.begin(), names.end(), meas.first);
                if(it != names.end()){
                    
                    
                    if (meas.first != pName){ // let's just deal with worms for now
                        for (auto & type : meas.second.object()){ //static vs dynamic
                            for (auto & func : type.second.object()){ //ij or ijkl of space
                                
                                auto const samples = type.first == "static" ?
                                    jInput["WL"][meas.first][func.first]["steps"].int64() :
                                    jInput["WL"][pName][pSpec]["steps"].int64();
                                
                                //need to do this ffor the other types
                                if (func.second.is<io::cvec>()){
                                    std::vector<ut::complex> data = jsx::at<io::cvec>(func.second);
                                    for (auto& x : data) x*=samples;
                                    jMeasurements(meas.first)(type.first)(func.first) << meas::fix(data, samples);
                                } else if (func.second.is<io::rvec>()) {
                                    std::vector<double> data =jsx::at<io::rvec>(func.second);
                                    for (auto& x : data) x*=samples;
                                    jMeasurements(meas.first)(type.first)(func.first) << meas::fix(data, samples);
                                }
                        
                            }
                        }
                    }
                        
                }
            }
         
            mpi::cout << "Done initializing measurements from previous run" << std::endl;
        }
    }
 
    
}




#endif
