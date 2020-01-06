#ifndef MEASUREMENTS
#define MEASUREMENTS

#include <vector>
#include <complex>
#include <fstream>
#include <valarray>
#include <cstring>
#include <random>

#include "../JsonX.h"
#include "../mpi/Utilities.h"
#include "../io/Vector.h"

//Achtung: es kann sein dass gewisse observabeln nicht gespeichert wurden, c.f. MonteCarlo.h

namespace meas {

    template<typename T>
    struct VecFix {
        inline static std::string name() { return name(T());};
        
        VecFix() = default;
        VecFix(VecFix const&) = default;
        VecFix(VecFix&&) = default;
        VecFix& operator=(VecFix const&) = default;
        VecFix& operator=(VecFix&&) = default;
        ~VecFix() = default;
        
        void add(std::vector<T> const& val, std::int64_t samples) {
            if(data_.size() == 0) data_.resize(val.size(), .0);
            if(data_.size() != val.size()) throw std::runtime_error(name() + "::add: missmatch in array size!");

            samples_ += samples; for(std::size_t n = 0; n < val.size(); ++n) data_[n] += val[n];
        };
        
        void write(jsx::value& dest) const {
            auto samples = samples_; mpi::reduce<mpi::op::sum>(samples, mpi::master);
            auto data = data_; mpi::reduce<mpi::op::sum>(data, mpi::master);
            
            if(mpi::rank() == mpi::master) {
                if(!samples) throw std::runtime_error(name() + "::write: no measurements taken !");  //Scheisse das sštt eh nie passiere.
                
                for(auto& x : data) x /= static_cast<double>(samples);
                
                data.b64() = true;  dest[io::Vector<T>::name()] = data;
            } else
                dest = "(rank " + std::to_string(mpi::rank()) + ") only the master collects the data";
        };
    private:
        std::int64_t samples_ = 0;
        io::Vector<T> data_;
        
        inline static std::string name(double const&) { return "meas::rvecfix";};
        inline static std::string name(std::complex<double> const&) { return "meas::cvecfix";};
    };


    template<typename T>
    struct VecVar {
        inline static std::string name() { return name(T());};
        
        VecVar() = default;
        VecVar(VecVar const&) = default;
        VecVar(VecVar&&) = default;
        VecVar& operator=(VecVar const&) = default;
        VecVar& operator=(VecVar&&) = default;
        ~VecVar() = default;
        
        void add(std::vector<T> const& val, std::int64_t samples) {
            if(val.size() > data_.size()) data_.resize(val.size(), .0);
            
            samples_ += samples; for(std::size_t n = 0; n < val.size(); ++n) data_[n] += val[n];
        };
    
        void write(jsx::value& dest) const {
            auto samples = samples_; mpi::reduce<mpi::op::sum>(samples, mpi::master);
            
            auto size = data_.size(); mpi::all_reduce<mpi::op::max>(size);
            auto data = data_; data.resize(size, .0); mpi::reduce<mpi::op::sum>(data, mpi::master);
            
            if(mpi::rank() == mpi::master) {
                if(!samples) throw std::runtime_error(name() + "::write: no measurements taken !");  //Scheisse das sštt eh nie passiere.
                
                for(auto& x : data) x /= static_cast<double>(samples);
                
                data.b64() = true; dest[io::Vector<T>::name()] = data;
            } else
                dest = "(rank " + std::to_string(mpi::rank()) + ") only the master collects the data";
        };
    private:
        std::int64_t samples_ = 0;
        io::Vector<T> data_;
        
        inline static std::string name(double const&) { return "meas::rvecvar";};
        inline static std::string name(std::complex<double> const&) { return "meas::cvecvar";};
    };

    
    template<typename T>
    struct FixSample {
        FixSample(T const& value, std::int64_t samples) : value(value), samples(samples) {};
        T const& value; std::int64_t samples;
    };
    
    template<typename T>
    struct VarSample {
        VarSample(T const& value, std::int64_t samples) : value(value), samples(samples) {};
        T const& value; std::int64_t samples;
    };
    
    
    template<typename T> inline FixSample<T> fix(T const& value, std::int64_t samples) { return FixSample<T>(value, samples);};
    template<typename T> inline VarSample<T> var(T const& value, std::int64_t samples) { return VarSample<T>(value, samples);};
    
    template<typename T>
    inline void operator<<(jsx::value& lhs, FixSample<T> const& rhs) {
        if(lhs.is<jsx::empty_t>()) lhs = VecFix<T>();
        lhs.at<VecFix<T>>().add(std::vector<T>(1, rhs.value), rhs.samples);
    }
    
    template<typename T>
    inline void operator<<(jsx::value& lhs, FixSample<std::vector<T>> const& rhs) {
        if(lhs.is<jsx::empty_t>()) lhs = VecFix<T>();
        lhs.at<VecFix<T>>().add(rhs.value, rhs.samples);
    }
    
    template<typename T>
    inline void operator<<(jsx::value& lhs, VarSample<std::vector<T>> const& rhs) {
        if(lhs.is<jsx::empty_t>()) lhs = VecVar<T>();
        lhs.at<VecVar<T>>().add(rhs.value, rhs.samples);
    }
    
    //--------------------------------------------------------------------------------------------------------------------------------
    
    inline void read_impl(jsx::value const& jMeasurements, double const sign, jsx::value& jObservables) {
        if(jMeasurements.is<jsx::object_t>())
            for(auto& entry : jMeasurements.object()) {
                if(entry.first == io::rvec::name()) {
                    auto vec = entry.second;
                    for(auto& x : jsx::at<io::rvec>(vec)) x /= sign;
                    jObservables = vec;
                } else if(entry.first == io::cvec::name()) {
                    auto vec = entry.second;
                    for(auto& x : jsx::at<io::cvec>(vec)) x /= sign;
                    jObservables = vec;
                } else
                    read_impl(entry.second, sign, jObservables[entry.first]);
            }
        
        if(jMeasurements.is<jsx::array_t>()) {
            jObservables = jsx::array_t(jMeasurements.size()); int index = 0;
            for(auto& elem : jMeasurements.array())
                read_impl(elem, sign, jObservables(index++));
        }
    };
    
    inline jsx::value read(jsx::value const& jMeasurements) {
        jsx::value jObservables;
        
        jsx::value jSign = jMeasurements("sign")(io::rvec::name());
        read_impl(jMeasurements, jsx::at<io::rvec>(jSign).at(0), jObservables);
        jObservables["sign"] = jsx::at<io::rvec>(jSign).at(0);
        
        return jObservables;
    };
}




#endif
