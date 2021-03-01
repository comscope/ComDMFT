#ifndef INCLUDE_MEASUREMENTS_ERROR_H
#define INCLUDE_MEASUREMENTS_ERROR_H

#include <vector>
#include <complex>
#include <cstring>

#include "../JsonX.h"
#include "../mpi/Utilities.h"
#include "../io/Vector.h"

//Achtung: es kann sein dass gewisse observabeln nicht gespeichert wurden, c.f. MonteCarlo.h

namespace meas {
    

    inline void error(std::vector<double>& arg, Jackknife) {
        auto temp = arg;  double const norm = mpi::number_of_workers();
        
        mpi::all_reduce<mpi::op::sum>(temp);
        for(std::size_t i = 0; i < arg.size(); ++i)
            temp[i] = (arg[i] - temp[i]/norm)*(arg[i] - temp[i]/norm);
        
        mpi::reduce<mpi::op::sum>(temp, mpi::master);
        for(std::size_t i = 0; i < arg.size(); ++i)
            arg[i] = 2.*std::sqrt((norm - 1)/norm*temp[i]);
    }
    
    
    template<typename E>
    inline void error(jsx::value& jArg, E) {
        if(jArg.is<io::rvec>())
            error(jArg.at<io::rvec>(), E());
        else if(jArg.is<io::cvec>()) {
            auto& arg = jArg.at<io::cvec>();
            
            std::vector<double> real, imag;
            for(auto& entry : arg) {
                real.push_back(entry.real()); imag.push_back(entry.imag());
            }
            error(real, E()); error(imag, E());
            for(std::size_t i = 0; i < arg.size(); ++i)
                arg[i] = {real[i], imag[i]};
        } else if(jArg.is<jsx::object_t>()) {
            for(auto& jEntry : jArg.object()) error(jEntry.second, E());
        } else if(jArg.is<jsx::array_t>()) {
            for(auto& jEntry : jArg.array()) error(jEntry, E());
        }
    }

    
    
    
    struct Error {
        Error() = default;
        ~Error() = default;
        
        void add(jsx::value const& jObservable, jsx::value const& jObservable0) {
            add(jMean_, jSquare_, jObservable, jObservable0);
        }
        jsx::value finalize(double norm, jsx::value const& jObservable0) {
            finalize(jMean_, jSquare_, jObservable0, norm);
            return std::move(jMean_);
        }
    private:
        jsx::value jMean_, jSquare_;
        
        static void add(jsx::value& jMean, jsx::value& jSquare, jsx::value const& jObservable, jsx::value const& jObservable0) {
            if(jObservable.is<io::rvec>()) {
                auto& obs = jObservable.at<io::rvec>();
                auto& obs0 = jObservable0.at<io::rvec>();
                
                if(!jMean.is<io::rvec>()) jMean = io::rvec(obs.size(), .0);
                if(!jSquare.is<io::rvec>()) jSquare = io::rvec(obs.size(), .0);
                
                auto& mean = jMean.at<io::rvec>();
                auto& square = jSquare.at<io::rvec>();
                
                for(std::size_t i = 0; i < obs.size(); ++i) {
                    mean[i] += obs[i]; square[i] += (obs[i] - obs0[i])*(obs[i] - obs0[i]);
                }
            } else if(jObservable.is<io::cvec>()) {
                auto& obs = jObservable.at<io::cvec>();
                auto& obs0 = jObservable0.at<io::cvec>();
                
                if(!jMean.is<io::cvec>()) jMean = io::cvec(obs.size(), .0);
                if(!jSquare.is<io::cvec>()) jSquare = io::cvec(obs.size(), .0);
                
                auto& mean = jMean.at<io::cvec>();
                auto& square = jSquare.at<io::cvec>();
                
                for(std::size_t i = 0; i < obs.size(); ++i) {
                    mean[i] += obs[i];
                    square[i] += std::complex<double>(
                        (obs[i].real() - obs0[i].real())*(obs[i].real() - obs0[i].real()),
                        (obs[i].imag() - obs0[i].imag())*(obs[i].imag() - obs0[i].imag())
                    );
                }
            } else if(jObservable.is<jsx::object_t>()) {
                for(auto& jEntry : jObservable.object())
                    add(jMean[jEntry.first], jSquare[jEntry.first], jEntry.second, jObservable0(jEntry.first));
            } else if(jObservable.is<jsx::array_t>()) {
                if(!jMean.is<jsx::array_t>()) jMean = jsx::array_t(jObservable.size());
                if(!jSquare.is<jsx::array_t>()) jSquare = jsx::array_t(jObservable.size());
                
                int index = 0;
                for(auto& jEntry : jObservable.array()) {
                    add(jMean[index], jSquare[index], jEntry, jObservable0(index)); ++index;
                }
            } else {
                jMean   = jObservable;
                jSquare = jObservable;
            }
        };
        
        static void finalize(jsx::value& jMean, jsx::value const& jSquare, jsx::value const& jObservable0, double const norm) {
            if(jMean.is<io::rvec>()) {
                auto& mean = jMean.at<io::rvec>();
                auto& square = jSquare.at<io::rvec>();
                auto& obs0 = jObservable0.at<io::rvec>();
                
                for(std::size_t i = 0; i < mean.size(); ++i) {
                    mean[i] = 2.*std::sqrt((norm - 1.)/norm*(square[i] - norm*(mean[i]/norm - obs0[i])*(mean[i]/norm - obs0[i])));
                }
            } else if(jMean.is<io::cvec>()) {
                auto& mean = jMean.at<io::cvec>();
                auto& square = jSquare.at<io::cvec>();
                auto& obs0 = jObservable0.at<io::cvec>();
                
                for(std::size_t i = 0; i < mean.size(); ++i)
                    mean[i] = {
                        2.*std::sqrt((norm - 1.)/norm*(square[i].real() - norm*(mean[i].real()/norm - obs0[i].real())*(mean[i].real()/norm - obs0[i].real()))),
                        2.*std::sqrt((norm - 1.)/norm*(square[i].imag() - norm*(mean[i].imag()/norm - obs0[i].imag())*(mean[i].imag()/norm - obs0[i].imag())))
                    };
            } else if(jMean.is<jsx::object_t>()) {
                for(auto& jEntry : jMean.object())
                    finalize(jEntry.second, jSquare(jEntry.first), jObservable0(jEntry.first), norm);
            } else if(jMean.is<jsx::array_t>()) {
                int index = 0;
                for(auto& jEntry : jMean.array()) {
                    finalize(jEntry, jSquare(index), jObservable0(index), norm); ++index;
                };
            }
        }
    };
}




#endif
