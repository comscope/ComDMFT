#ifndef CTQMC_INCLUDE_DATA_H
#define CTQMC_INCLUDE_DATA_H

#include <vector>
#include <tuple>

#include "Utilities.h"
#include "impurity/Diagonal.h"
#include "impurity/Operators.h"
#include "impurity/Dynamic.h"
#include "impurity/Tensor.h"
#include "bath/Hyb.h"

#include "impurity/Observables.h"

#include "config/Worms.h"

#include "../../include/mpi/Utilities.h"
#include "../../include/atomic/Generate.h"
#include "../../include/options/Options.h"

namespace data {
    
    namespace impl {
        
        template<typename... Types> struct get_tuple_type_index;
        
        template<typename T, typename... Types>
        struct get_tuple_type_index<T, std::tuple<Types...>> {
            constexpr static std::size_t value = ut::get_index_by_type<T, Types...>::value;
        };
        
    }
    
    
    template<typename Value>
    using Optional = std::tuple<
    std::unique_ptr<imp::itf::BullaOperators<Value>>,
    std::unique_ptr<imp::itf::Occupation<Value>>,
    std::unique_ptr<imp::itf::BullaOccupation<Value>>,
    std::unique_ptr<imp::Tensor<Value>>
    >;
    
    
    template<typename Value>
    struct Data {
        Data() = delete;
        template<typename Mode>
        Data(jsx::value const& jParams, Mode) {
            ut::beta = ut::Beta(jParams("beta").real64());         // Bad Bad Bad ...

            filling_ = jsx::at<io::rvec>(jParams("hloc")("filling")); filling_.insert(filling_.begin(), 0);
            if(jParams.is("dyn")) dyn_.reset(new imp::Simple(jParams, jParams("dyn")));
            eig_.reset(new imp::EigenValues<Mode>(jParams, jParams("hloc")("eigen values"), filling(), dyn()));
            ide_.reset(new imp::Operator<Mode, Value>('1', eig()));
            ops_.reset(new imp::Operators<Mode, Value>(jParams, jParams("operators"), eig()));
            hyb_.reset(new bath::Hyb<Value>(jParams, jParams("hybridisation")("matrix"), jParams("hybridisation")("functions")));
        }
        Data(Data const&) = delete;
        Data(Data&&) = delete;
        Data& operator=(Data const&) = delete;
        Data& operator=(Data&&) = delete;
        ~Data() = default;

        std::vector<double> const& filling() const {
            return filling_;
        };
        imp::itf::EigenValues const& eig() const {
            return *eig_;
        };
        imp::itf::Operator<Value> const& ide() const {
            return *ide_;
        };
        imp::itf::Operators<Value> const& ops() const {
            return *ops_;
        };
        bath::Hyb<Value> const& hyb() const {
            return *hyb_;
        };
        imp::Simple const* dyn() const {
            return dyn_.get();
        };

        
        template<typename T> std::unique_ptr<T>& opt() {
            return std::get<impl::get_tuple_type_index<std::unique_ptr<T>, Optional<Value>>::value>(optional_);
        };
        template<typename T> T const& opt() const {
            return *std::get<impl::get_tuple_type_index<std::unique_ptr<T>, Optional<Value>>::value>(optional_).get();
        };

    private:
        std::vector<double> filling_;
        std::unique_ptr<imp::itf::EigenValues> eig_;
        std::unique_ptr<imp::itf::Operator<Value>> ide_;
        std::unique_ptr<imp::itf::Operators<Value>> ops_;
        std::unique_ptr<bath::Hyb<Value>> hyb_;
        std::unique_ptr<imp::Simple> dyn_;
        
        Optional<Value> optional_;
    };
    
    
    template<typename Mode, typename Value>
    struct init_worm_data_functor {
        template<typename W>
        void operator()(ut::wrap<W> w, jsx::value const& jParams, data::Data<Value>& data) const {
            if(jParams.is(W::name()))
                cfg::init_worm_data<Mode>(w, jParams, data);
        }
    };
    
    template<typename Mode, typename Value>
    void setup_data(jsx::value const& jParams, data::Data<Value>& data)
    {
        cfg::for_each_type<cfg::Worm>::apply(init_worm_data_functor<Mode, Value>(), jParams, data);
    };
    
}

#endif
