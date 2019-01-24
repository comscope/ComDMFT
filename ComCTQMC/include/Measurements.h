#ifndef MEASUREMENTS
#define MEASUREMENTS

#include <vector>
#include <complex>
#include <fstream>
#include <valarray>
#include <cstring>
#include <random>

#include "MPIUtilities.h"
#include "JsonX.h"
#include "IO.h"

//Achtung: es kann sein dass gewisse observabeln nicht gespeichert wurden, c.f. MonteCarlo.h

namespace meas {

    template<typename T, typename std::enable_if<std::is_same<double, T>::value || std::is_same<std::complex<double>, T>::value, int>::type = 0>
    struct VectorObs {
    protected:
        std::int64_t samples_ = 0;
        io::Vector<T> data_;
        
        VectorObs() = default;
        VectorObs(VectorObs const&) = default;
        VectorObs(VectorObs&&) = default;
        VectorObs& operator=(VectorObs const&) = default;
        VectorObs& operator=(VectorObs&&) = default;
        
        void add(std::vector<T> const& val, std::int64_t samples) {
            if(data_.size() != val.size())
                throw std::runtime_error("meas::VectorObs::add: missmatch in array size!");
            
            samples_ += samples;
            for(std::size_t n = 0; n < data_.size(); ++n) data_[n] += val[n];
        };
        
        static void write(VectorObs<T> const& source, jsx::value& dest) {
            VectorObs<T> reduced;
            
#ifdef HAVE_MPI
            MPI_Reduce(&source.samples_, mpi::rank() == mpi::master ? &reduced.samples_ : 0, 1, MPI_INT64_T, MPI_SUM, mpi::master, MPI_COMM_WORLD);
            
            if(mpi::rank() == mpi::master) reduced.data_.resize(source.data_.size());

            MPI_Reduce(source.data_.data(), mpi::rank() == mpi::master ? reduced.data_.data() : 0, source.data_.size(), std::is_same<double, T>::value ? MPI_DOUBLE : MPI_C_DOUBLE_COMPLEX, MPI_SUM, mpi::master, MPI_COMM_WORLD);
#else
            reduced = source;
#endif
            
            if(mpi::rank() == mpi::master) {
                if(!reduced.samples_)
                    throw std::runtime_error("meas::VectorObs::write: no measurements taken !");  //Scheisse das sštt eh nie passiere.
                
                for(auto& x : reduced.data_) x /= static_cast<double>(reduced.samples_);
                
                reduced.data_.write(dest[std::is_same<double, T>::value ? "__rmeas__" : "__cmeas__"]);
            } else
                dest = "(rank " + std::to_string(mpi::rank()) + ") only the master collects the data";
        };
    };
    
    
    template<typename T>
    struct VecObsFix : VectorObs<T> {
        VecObsFix() = default;
        VecObsFix(VecObsFix const&) = default;
        VecObsFix(VecObsFix&&) = default;
        VecObsFix& operator=(VecObsFix const&) = default;
        VecObsFix& operator=(VecObsFix&&) = default;
        
        void add(std::vector<T> const& val, std::int64_t samples) {
            if(this->data_.size() == 0) this->data_.resize(val.size(), .0);
            VectorObs<T>::add(val, samples);
        };
        
        void write(jsx::value& arg) const {
            VectorObs<T>::write(*this, arg);
        };
    };
    
    
    template<typename T>
    struct VecObsVar : VectorObs<T> {
        VecObsVar() = default;
        VecObsVar(VecObsVar const&) = default;
        VecObsVar(VecObsVar&&) = default;
        VecObsVar& operator=(VecObsVar const&) = default;
        VecObsVar& operator=(VecObsVar&&) = default;
        
        void add(std::vector<T> const& val, std::int64_t samples) {
            if(val.size() > this->data_.size()) this->data_.resize(val.size(), .0);
            VectorObs<T>::add(val, samples);
        };
        
        void write(jsx::value& arg) const {
#ifdef HAVE_MPI
            std::int64_t this_size = this->data_.size(), max_size;
            MPI_Allreduce(&this_size, &max_size, 1, MPI_INT64_T, MPI_MAX, MPI_COMM_WORLD);
            VecObsVar<T> resized = *this; resized.data_.resize(max_size, .0);
            
            VectorObs<T>::write(resized, arg);
#else
            VectorObs<T>::write(*this, arg); 
#endif
        };
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
    
    
    template<typename T> FixSample<T> fix(T const& value, std::int64_t samples) { return FixSample<T>(value, samples);};
    template<typename T> VarSample<T> var(T const& value, std::int64_t samples) { return VarSample<T>(value, samples);};
    
    //--------------------------------------------------------------------------------------------------------------------------------
    
    inline void read_impl(jsx::value const& jMeasurements, double const sign, jsx::value& jObservables) {
        if(jMeasurements.type() == jsx::type::object)
            for(auto& entry : jMeasurements.object()) {
                if(entry.first == "__rmeas__") {
                    auto vec = entry.second;
                    for(auto& x : jsx::at<io::rvec>(vec)) x /= sign;
                    jObservables = vec;
                } else if(entry.first == "__cmeas__") {
                    auto vec = entry.second;
                    for(auto& x : jsx::at<io::cvec>(vec)) x /= sign;
                    jObservables = vec;
                } else
                    read_impl(entry.second, sign, jObservables[entry.first]);
            }
        
        if(jMeasurements.type() == jsx::type::array) {
            jObservables = jsx::array(jMeasurements.size()); int index = 0;
            for(auto& elem : jMeasurements.array())
                read_impl(elem, sign, jObservables(index++));
        }
    };
    
    inline jsx::value read(jsx::value const& jMeasurements) {
        jsx::value jObservables;
        
        jsx::value jSign = jMeasurements("Sign")("__rmeas__");
        read_impl(jMeasurements, jsx::at<io::rvec>(jSign).at(0), jObservables);
        jObservables["Sign"] = jsx::at<io::rvec>(jSign).at(0);
        
        return jObservables;
    };
}


template<typename T>
void operator<<(jsx::value& lhs, meas::FixSample<T> const& rhs) {
    if(lhs.type() == jsx::type::empty) lhs = meas::VecObsFix<T>();
    lhs.at<meas::VecObsFix<T>>().add(std::vector<T>(1, rhs.value), rhs.samples);
}

template<typename T>
void operator<<(jsx::value& lhs, meas::FixSample<std::vector<T>> const& rhs) {
    if(lhs.type() == jsx::type::empty) lhs = meas::VecObsFix<T>();
    lhs.at<meas::VecObsFix<T>>().add(rhs.value, rhs.samples);
}

template<typename T>
void operator<<(jsx::value& lhs, meas::VarSample<std::vector<T>> const& rhs) {
    if(lhs.type() == jsx::type::empty) lhs = meas::VecObsVar<T>();
    lhs.at<meas::VecObsVar<T>>().add(rhs.value, rhs.samples);
}

#endif
