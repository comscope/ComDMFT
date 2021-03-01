#ifndef INCLUDE_MPI_BASIC_H
#define INCLUDE_MPI_BASIC_H

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <string>
#include <vector>
#include <complex>

//--------------------------------------------------schö tö tiä tü mö tiä par la barbischätöööötötötötö-------------------------------------------------------

namespace mpi {
    
    enum : unsigned {
        empty       = 0x0,
        fundamental = 0x1,
        arithmetic  = 0x2,
        ordered     = 0x4,
        logical     = 0x8,
        bitwise     = 0x10
    };
    
    
    template<unsigned p>
    struct has_property {
        enum : unsigned { value = p };
    };
    
    
    template<typename T> struct property              :  has_property< empty > {};

    template<> struct property<char>                  :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<signed char>           :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<unsigned char>         :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<short>                 :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<unsigned short>        :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<int>                   :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<unsigned>              :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<long>                  :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<unsigned long>         :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<long long>             :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<unsigned long long>    :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<float>                 :  has_property< fundamental | arithmetic | ordered > {};
    template<> struct property<double>                :  has_property< fundamental | arithmetic | ordered > {};
    template<> struct property<std::complex<float>>   :  has_property< fundamental | arithmetic > {};
    template<> struct property<std::complex<double>>  :  has_property< fundamental | arithmetic > {};
    
    
    template<typename T, unsigned m>
    struct data_compatible_if {
        enum : bool { value = ( (property<T>::value & m) == m ) };
    };
    
    
    namespace op {
        struct min {};  struct max {}; struct sum {};  struct prod {};
        struct land {}; struct lor {}; struct lxor {};
        struct band {}; struct bor {}; struct bxor {};
    };
    
    
    template<typename T, typename Op > struct data_op_compatible { enum : bool { value = false }; };
    
    template<typename T> struct data_op_compatible<T, op::min>   :  data_compatible_if< T, fundamental | ordered > {};
    template<typename T> struct data_op_compatible<T, op::max>   :  data_compatible_if< T, fundamental | ordered > {};
    template<typename T> struct data_op_compatible<T, op::sum>   :  data_compatible_if< T, fundamental | arithmetic > {};
    template<typename T> struct data_op_compatible<T, op::prod>  :  data_compatible_if< T, fundamental | arithmetic > {};
    template<typename T> struct data_op_compatible<T, op::land>  :  data_compatible_if< T, fundamental | logical > {};
    template<typename T> struct data_op_compatible<T, op::lor>   :  data_compatible_if< T, fundamental | logical > {};
    template<typename T> struct data_op_compatible<T, op::lxor>  :  data_compatible_if< T, fundamental | logical > {};
    template<typename T> struct data_op_compatible<T, op::band>  :  data_compatible_if< T, fundamental | bitwise > {};
    template<typename T> struct data_op_compatible<T, op::bor>   :  data_compatible_if< T, fundamental | bitwise > {};
    template<typename T> struct data_op_compatible<T, op::bxor>  :  data_compatible_if< T, fundamental | bitwise > {};
    
    
#ifdef HAVE_MPI
    MPI_Op get_op(op::min const&)   { return MPI_MIN; };
    MPI_Op get_op(op::max const&)   { return MPI_MAX; };
    MPI_Op get_op(op::sum const&)   { return MPI_SUM; };
    MPI_Op get_op(op::prod const&)  { return MPI_PROD; };
    MPI_Op get_op(op::land const&)  { return MPI_LAND; };
    MPI_Op get_op(op::lor const&)   { return MPI_LOR; };
    MPI_Op get_op(op::lxor const&)  { return MPI_LXOR; };
    MPI_Op get_op(op::band const&)  { return MPI_BAND; };
    MPI_Op get_op(op::bor const&)   { return MPI_BOR; };
    MPI_Op get_op(op::bxor const&)  { return MPI_BXOR; };
    
    template<typename T> MPI_Datatype get_data_type(T const&); // redundant because of sfinae, however, sometimes redundancy is good ...
    
    MPI_Datatype get_data_type(char const&)                 { return MPI_CHAR; };
    MPI_Datatype get_data_type(signed char const&)          { return MPI_SIGNED_CHAR; };
    MPI_Datatype get_data_type(unsigned char const&)        { return MPI_UNSIGNED_CHAR; };
    MPI_Datatype get_data_type(short const&)                { return MPI_SHORT; };
    MPI_Datatype get_data_type(unsigned short const&)       { return MPI_UNSIGNED_SHORT; };
    MPI_Datatype get_data_type(int const&)                  { return MPI_INT; };
    MPI_Datatype get_data_type(unsigned const&)             { return MPI_UNSIGNED; };
    MPI_Datatype get_data_type(long const&)                 { return MPI_LONG; };
    MPI_Datatype get_data_type(unsigned long const&)        { return MPI_UNSIGNED_LONG; };
    MPI_Datatype get_data_type(long long const&)            { return MPI_LONG_LONG; };
    MPI_Datatype get_data_type(unsigned long long const&)   { return MPI_UNSIGNED_LONG_LONG; };
    MPI_Datatype get_data_type(float const&)                { return MPI_FLOAT; };
    MPI_Datatype get_data_type(double const&)               { return MPI_DOUBLE; };
    MPI_Datatype get_data_type(std::complex<float> const&)  { return MPI_C_FLOAT_COMPLEX; };
    MPI_Datatype get_data_type(std::complex<double> const&) { return MPI_C_DOUBLE_COMPLEX; };
#endif
    
    
    int const master = 0;
    
    inline int rank() {
        int temp = 0;
#ifdef HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &temp);
#endif
        return temp;
    };
    
    inline int number_of_workers() {
        int temp = 1;
#ifdef HAVE_MPI
        MPI_Comm_size(MPI_COMM_WORLD, &temp);
#endif
        return temp;
    };
    
    inline int processor_name_size() {
#ifdef HAVE_MPI
        return MPI_MAX_PROCESSOR_NAME;
#else
        return 1;
#endif
    }

    inline std::vector<char> processor_name() {
#ifdef HAVE_MPI 
        std::vector<char> name(processor_name_size(), '\0'); int length;
        MPI_Get_processor_name(name.data(), &length);
        return name;
#else
        return { '\0' };
#endif
    }
    
    inline void barrier() {
#ifdef HAVE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    };
    
    
    template<typename Op, typename T, typename std::enable_if<data_op_compatible<T, Op>::value, int>::type = 0>
    void reduce(T& arg, int root) {
#ifdef HAVE_MPI
        T result; MPI_Reduce(&arg, &result, 1, get_data_type(T()), get_op(Op()), root, MPI_COMM_WORLD);
        if(rank() == root) arg = result;
#endif
    };
    
    template<typename Op, typename T, typename std::enable_if<data_op_compatible<T, Op>::value, int>::type = 0>
    void reduce(std::vector<T>& arg, int root) {
#ifdef HAVE_MPI
        std::vector<T> result(rank() == root ? arg.size() : 0);
        MPI_Reduce(arg.data(), result.data(), arg.size(), get_data_type(T()), get_op(Op()), root, MPI_COMM_WORLD);
        if(rank() == root) arg = std::move(result);
#endif
    };
    
    
    template<typename Op, typename T, typename std::enable_if<data_op_compatible<T, Op>::value, int>::type = 0>
    void all_reduce(T& arg) {
#ifdef HAVE_MPI
        T result; MPI_Allreduce(&arg, &result, 1, get_data_type(T()), get_op(Op()), MPI_COMM_WORLD);
        arg = result;
#endif
    };
    
    template<typename Op, typename T, typename std::enable_if<data_op_compatible<T, Op>::value, int>::type = 0>
    void all_reduce(std::vector<T>& arg) {
#ifdef HAVE_MPI
        std::vector<T> result(arg.size());
        MPI_Allreduce(arg.data(), result.data(), arg.size(), get_data_type(T()), get_op(Op()), MPI_COMM_WORLD);
        arg = std::move(result);
#endif
    };
 
    
    template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type = 0>
    void bcast(T& arg, int root) {
#ifdef HAVE_MPI
        MPI_Bcast(&arg, 1, get_data_type(T()), root, MPI_COMM_WORLD);
#endif
    };
    
    template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type = 0>
    void bcast(std::vector<T>& arg, int root) {
#ifdef HAVE_MPI
        MPI_Bcast(arg.data(), arg.size(), get_data_type(T()), root, MPI_COMM_WORLD);
#endif
    };
    
    template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type = 0>
    void bcast(std::basic_string<T>& arg, int root) {
#ifdef HAVE_MPI
        MPI_Bcast(&arg.front(), arg.size(), get_data_type(T()), root, MPI_COMM_WORLD);
#endif
    };
    
    
    template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type = 0>
    void gather(std::vector<T>& arg, int root) {
#ifdef HAVE_MPI
        std::vector<T> result(rank() == root ? arg.size()*mpi::number_of_workers() : 0);
        MPI_Gather(arg.data(), arg.size(), get_data_type(T()), result.data(), arg.size(), get_data_type(T()), root, MPI_COMM_WORLD);
        if(rank() == root) arg = std::move(result);
#endif
    };
    
    
    template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type = 0>
    void scatter(std::vector<T>& arg, int size, int root) {
        if(rank() == root && size*number_of_workers() != static_cast<int>(arg.size()))
            throw std::runtime_error("mpi::scatter: send buffer has wrong size");        
#ifdef HAVE_MPI
        std::vector<T> result(size);
        MPI_Scatter(arg.data(), size, get_data_type(T()), result.data(), size, get_data_type(T()), root, MPI_COMM_WORLD);
        arg = std::move(result);
#endif
    };
    
}


#endif
