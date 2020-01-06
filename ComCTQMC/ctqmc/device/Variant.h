#ifndef VARIANT
#define VARIANT

// is there a way around ugly host - device namespace without getting warnings from compiler ?

namespace var {
    
    template<typename... Args> union variant {};
    
    template<typename T, typename... Args>
    union variant<T, Args...> {
        T value;
        variant<Args...> next;
    };
    

    template<typename T, typename U, typename... Args>
    struct get_impl {
        __host__ __device__ static T& at(variant<U, Args...>& arg) { return get_impl<T, Args...>::at(arg.next);};
    };
    
    template<typename T, typename... Args>
    struct get_impl<T, T, Args...> {
        __host__ __device__ static T& at(variant<T, Args...>& arg) { return arg.value;};
    };
    
    namespace host {
        template<typename T, typename... Args>
        T& get(variant<Args...>& arg) {
            return get_impl<T, Args...>::at(arg);
        };
    }

    namespace device {
        template<typename T, typename... Args>
        __device__ T& get(variant<Args...>& arg) {
            return get_impl<T, Args...>::at(arg);
        };
    }
    
    
    template<typename T, typename V> struct index;
    
    template<typename T, typename U, typename... Args>
    struct index<T, variant<U, Args...>> {
        static constexpr int value = index<T, variant<Args...>>::value + 1;
    };
    
    template<typename T, typename... Args>
    struct index<T, variant<T, Args...>> {
        static constexpr int value = 0;
    };
    
}


#endif
