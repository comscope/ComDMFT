#ifndef INCLUDE_JSONX
#define INCLUDE_JSONX

#include <iostream>
#include <cstring>
#include <vector>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <cstdlib>
#include <cerrno>
#include <utility>
#include <limits>
#include <map>


namespace jsx {
    
    struct value;
    
    template<typename T> struct trait {
        constexpr static bool is_json = false;
        static std::string name() {
            return T::name();
        };
        static void to_json(T const& t, value& dest) {
            t.write(dest);
        };
    };
    
    struct empty {};  struct null {};
    
// ugly macro nonsense

#define COMMA ,
    
#define JSON_TRAIT(TYPE, JSX_TYPE )                            \
using JSX_TYPE = TYPE;                                         \
template<> struct trait<JSX_TYPE> {                            \
    constexpr static bool is_json = true;                      \
    static std::string name() {                                \
        return # JSX_TYPE;                                     \
    };                                                         \
    static void to_json(JSX_TYPE const& t, value& dest) {      \
        throw std::runtime_error("jsx: is already json type"); \
    };                                                         \
};                                                             
    
    JSON_TRAIT(empty,                             empty_t  );
    JSON_TRAIT(null,                              null_t   );
    JSON_TRAIT(bool,                              boolean_t);
    JSON_TRAIT(std::int64_t,                      int64_t  );
    JSON_TRAIT(double,                            real64_t );
    JSON_TRAIT(std::string,                       string_t );
    JSON_TRAIT(std::vector<value>,                array_t  );
    JSON_TRAIT(std::map<std::string COMMA value>, object_t );
    
#undef COMMA
#undef JSON_TRAIT

// otherwise fancy initilasation does not work (whatever the name of this nonsense is, initialiser list maybe ???) 

    template<typename T> struct map_trait    { using type = T; };
    
    template<> struct map_trait<char const*> { using type = string_t; };
    template<> struct map_trait<int>         { using type = int64_t; };
    
    
    struct value {
        value() {
            manage_t<empty_t>::create(data_);
        };
        value(value const& other) {
            other.data_.manage(Op::clone, other.data_, &this->data_);
        };
        value(value&& other) noexcept {
            data_ = other.data_; manage_t<empty_t>::create(other.data_);
        };
        value& operator=(value const& other) {            //self-assignement safe
            Data temp; other.data_.manage(Op::clone, other.data_, &temp);
            data_.manage(Op::destroy, data_, nullptr); data_ = temp; return *this;
        };
        value& operator=(value&& other) noexcept {        //self-assignement safe
            Data temp = other.data_; manage_t<empty_t>::create(other.data_);
            data_.manage(Op::destroy, data_, nullptr); data_ = temp; return *this;
        };
        ~value() {
            data_.manage(Op::destroy, data_, nullptr);
        };
        
        template<typename T, typename Td = typename map_trait<typename std::decay<T>::type>::type, typename = typename std::enable_if<!std::is_same<value, Td>::value>::type>
        value(T&& t) {
            manage_t<Td>::create(data_, std::forward<T>(t));
        }
        
        template<typename T, typename Td = typename map_trait<typename std::decay<T>::type>::type, typename = typename std::enable_if<!std::is_same<value, Td>::value>::type>
        value& operator=(T&& t) {
            data_.manage(Op::destroy, data_, nullptr);  manage_t<Td>::create(data_, std::forward<T>(t)); return *this;
        }

        
        template<typename T>
        bool is() const {
            return data_.manage == &manage_t<T>::apply;
        }
        
        template<typename T> T& at() {
            return at_impl<T>();
        }
        template<typename T> T const& at() const {
            return at_impl<T>();
        }
    
        
        boolean_t& boolean() { return at<boolean_t>(); };
        string_t&  string()  { return at<string_t>(); };
        array_t&   array()   { return at<array_t>(); };
        object_t&  object()  { return at<object_t>(); };
        
        int64_t int64() const {
            if(is<real64_t>() && static_cast<real64_t>(static_cast<int64_t>(unsafe_at_impl<real64_t>())) == unsafe_at_impl<real64_t>()) return unsafe_at_impl<real64_t>();
            return at<int64_t>();
        };
        real64_t real64() const {
            if(is<int64_t>() && static_cast<int64_t>(static_cast<real64_t>(unsafe_at_impl<int64_t>())) == unsafe_at_impl<int64_t>()) return unsafe_at_impl<int64_t>();
            return at<real64_t>();
        };
        boolean_t const& boolean() const { return at<boolean_t>(); };
        string_t  const& string()  const { return at<string_t>(); };
        array_t   const& array()   const { return at<array_t>(); };
        object_t  const& object()  const { return at<object_t>(); };
    
        
        std::size_t size() const { //eleganz vo arsch vo chamel aber praktisch
            if(is<array_t>()) return array().size();
            if(is<object_t>()) return object().size();
            
            throw std::runtime_error("jsx::value: found " + name() + " type instead of jsx::array_t or jsx::object_ type for size request.");
        };
        
        value& operator[](std::size_t index) {
            return array()[index];
        };
        value& operator()(std::size_t index) {
            if(!(index < array().size())) throw std::runtime_error("jsx::value: invalid array index.");
            return array()[index];
        };
        value const& operator()(std::size_t index) const {
            if(!(index < array().size())) throw std::runtime_error("jsx::value: invalid array index.");
            return array()[index];
        };
        
        std::size_t is(std::string const& key) const {
            return object().count(key);
        };
        value& operator[](std::string const& key) {
            if(!is<object_t>()) *this = object_t();
            return object()[key];
        };
        value& operator()(std::string const& key) {
            if(!is(key)) throw std::runtime_error("jsx::value: key \"" + key + "\" not found.");
            return object()[key];
        };
        value const& operator()(std::string const& key) const {
            if(!is(key)) throw std::runtime_error("jsx::value: key \"" + key + "\" not found.");
            return object().at(key);
        };
        
        bool is_json() const {
            bool is; data_.manage(Op::is_json, data_, &is); return is;
        };
        value to_json() const {
            value v; data_.manage(value::Op::to_json, data_, &v); return v;
        };
        std::string name() const {
            std::string name; data_.manage(Op::name, data_, &name); return name;
        };
        
    private:
        enum class Op { clone, destroy, is_json, to_json, name };
        
        struct Data {
            union {
                void* ptr;
                std::aligned_storage<sizeof(void*), alignof(void*)>::type raw;
            };
            void (*manage) (Op, Data const&, void*);
        } data_;
        
        
        template<typename T>
        struct manage_ptr {
            template<typename... Args>
            static void create(Data& data, Args&&... args) {
                data.ptr = new T(std::forward<Args>(args)...);
                data.manage = &manage_ptr<T>::apply;
            }
            static void* get_mem(Data const& data) {
                return data.ptr;
            };
            static void apply(Op op, Data const& data, void* arg) {
                switch(op) {
                    case Op::clone:
                        static_cast<Data*>(arg)->ptr = new T(*static_cast<T const*>(get_mem(data)));
                        static_cast<Data*>(arg)->manage = &manage_ptr<T>::apply; break;
                    case Op::destroy:
                        delete static_cast<T*>(get_mem(data)); break;
                    case Op::is_json:
                        *static_cast<bool*>(arg) = trait<T>::is_json; break;
                    case Op::to_json:
                        trait<T>::to_json(*static_cast<T const*>(get_mem(data)), *static_cast<value*>(arg)); break;
                    case Op::name:
                        *static_cast<std::string*>(arg) = trait<T>::name(); break;
                }
            };
        };
        
        template<typename T>
        struct manage_raw {
            template<typename... Args>
            static void create(Data& data, Args&&... args) {
                new(&data.raw) T(std::forward<Args>(args)...);
                data.manage = &manage_raw<T>::apply;
            }
            static void* get_mem(Data& data) {
                return &data.raw;
            };
            static void const* get_mem(Data const& data) {
                return &data.raw;
            };
            static void apply(Op op, Data const& data, void* arg) {
                switch(op) {
                    case Op::clone:
                        new(&static_cast<Data*>(arg)->raw) T(*static_cast<T const*>(get_mem(data)));
                        static_cast<Data*>(arg)->manage = &manage_raw<T>::apply; break;
                    case Op::destroy:
                        static_cast<T const*>(get_mem(data))->~T(); break;
                    case Op::is_json:
                        *static_cast<bool*>(arg) = trait<T>::is_json; break;
                    case Op::to_json:
                        trait<T>::to_json(*static_cast<T const*>(get_mem(data)), *static_cast<value*>(arg)); break;
                    case Op::name:
                        *static_cast<std::string*>(arg) = trait<T>::name(); break;
                }
            };
        };
        
        template<typename T> using manage_t = typename std::conditional<sizeof(T) <= sizeof(void*) && alignof(void*)%alignof(T) == 0, manage_raw<T>, manage_ptr<T>>::type;
        
        
        template<typename T> T& unsafe_at_impl() {
            return *static_cast<T*>(manage_t<T>::get_mem(data_));
        }
        template<typename T> T const& unsafe_at_impl() const {
            return *static_cast<T const*>(manage_t<T>::get_mem(data_));
        }
        
        template<typename T>
        T& at_impl() {
            if(data_.manage == &manage_t<T>::apply) return unsafe_at_impl<T>();
            throw std::runtime_error("found " + name() + " type instead of " + trait<T>::name() + " type");
        }
        template<typename T>
        T const& at_impl() const {
            if(data_.manage == &manage_t<T>::apply) return unsafe_at_impl<T>();
            throw std::runtime_error("found " + name() + " type instead of " + trait<T>::name() + " type");
        }
        
    };
    
    
    template<typename T>
    T& at(value& source) {
        if(source.is_json()) {
            T dest; dest.read(source); source = std::move(dest);
        }
        return source.at<T>();
    };
    template<typename T>
    T const& at(value const& source) {
        return source.at<T>();
    };
    
    
    
    
    inline char const* parse(char const*, value&);
    
    inline const char* skip_white_space(const char* it) {
        while(*it != '\0' && std::isspace(*it)) ++it;
        return it;
    }
    
    inline char const* parse(char const* it, std::string& s) {
        ++it;
        
        char const* end = it;
        while(*end != '\0' && *end != '\"') ++end;
        if(*end == '\0' ) throw std::runtime_error("jsx::parser: error while reading string");
        s.append(it, end);
        if(*(end - 1) != '\\') return ++end;
        
        s.append(end, end + 1);      // keep fucking escaped double quotes ...
        return parse(end, s);
    }
    
    inline char const* parse(char const* it, array_t& a) {
        it = skip_white_space(++it);
        if(*it == ']') return ++it;
        
        while(1) {
            a.push_back(value());
            it = parse(it, a.back());
            if(*it != ',') break;
            ++it;
        };
        
        if(*it != ']') throw std::runtime_error("jsx::parser: error while reading array");
        return ++it;
    }
    
    inline char const* parse(char const* it, object_t& o) {
        it = skip_white_space(++it);
        if(*it == '}') return ++it;
        
        while(1) {
            if(*it != '\"') throw std::runtime_error("jsx::parser: object key not valid");
            std::string key; it = parse(it, key);
            
            it = skip_white_space(it);
            if(*it != ':') throw std::runtime_error("jsx::parser: no object value found");
            it = parse(++it, o[key]);
            
            if(*it != ',') break;
            it = skip_white_space(++it);
        };
        
        if(*it != '}') throw std::runtime_error("jsx::parser: error while reading object");
        return ++it;
    }
    
    constexpr long long int lli_int64_max = std::numeric_limits<int64_t>::max();
    constexpr long long int lli_int64_min = std::numeric_limits<int64_t>::min();
    
    inline char const* parse(char const* it, value& v) {
        it = skip_white_space(it);
        
        switch(*it) {
            case '{':
                v = object_t();
                it = parse(it, v.object());
                break;
            case '[':
                v = array_t();
                it = parse(it, v.array());
                break;
            case '\"':
                v = string_t();
                it = parse(it, v.string());
                break;
            case 'n':
                if(*++it == 'u' && *++it == 'l' && *++it == 'l') v = null_t(); else throw std::runtime_error("jsx::parser: error while reading null");
                ++it; break;
            case 't':
                if(*++it == 'r' && *++it == 'u' && *++it == 'e') v = true; else throw std::runtime_error("jsx::parser: error while reading true");
                ++it; break;
            case 'f':
                if(*++it == 'a' && *++it == 'l' && *++it == 's' && *++it == 'e') v = false; else throw std::runtime_error("jsx::parser: error while reading true");
                ++it; break;
            default:
                errno = 0; char* lli_end; long long int lli = strtoll(it, &lli_end, 10);
                int int64_fail = errno != 0 || lli_end == it || lli_int64_max < lli || lli < lli_int64_min;
                
                errno = 0; char* d_end; double d = strtod(it, &d_end);
                int d_fail = errno != 0 || d_end == it;
                
                if(!int64_fail && lli_end == d_end) {
                    v = static_cast<int64_t>(lli); it = lli_end; break;
                }
                
                if(!d_fail) {
                    v = d; it = d_end; break;
                }
                
                throw std::runtime_error("jsx::parser: empty found !");
        }
        
        return skip_white_space(it);
    }
    
    
    inline value read(std::string const name) {
        std::ifstream file(name.c_str());  std::string buffer;
        
        if(!file) throw std::runtime_error("mpi: file " + name + " not found !");
        
        file.seekg(0, std::ios::end);
        int size = file.tellg();
        file.seekg(0, std::ios::beg);
        
        buffer.resize(size + 1);
        file.read(&buffer.front(), size);
        buffer[size] = '\0';
        
        file.close();

        value value;
        if(*parse(&buffer.front(), value) != '\0') throw std::runtime_error("jsx::parser: stream contains more ...");
        return value;
    }
    
    inline value read(std::string const name, value dvalue) {
        if(!std::ifstream(name.c_str())) return dvalue;
        return read(name);
    }
    
    
    inline void write(value const&, std::ostream&, int, int, bool);
    
    inline void write(array_t const& a, std::ostream& stream, int indent = 4, int pos = 0, bool in_array = false) {
        if(in_array) stream << '\n' << std::string(pos += indent, ' ');
        stream << "[ ";
        auto it = a.begin();
        if(it != a.end())
            while(1) {
                write(*it, stream, indent, pos, true);
                if(++it != a.end())
                    stream << ", ";
                else
                    break;
            }
        stream << " ]";
    }
    
    inline void write(object_t const& o, std::ostream& stream, int indent = 4, int pos = 0, bool in_array = false) {
        stream << '{';
        auto it = o.begin();
        if(it != o.end())
            while(1) {
                stream << '\n' << std::string(pos + indent, ' ') << '\"' << it->first << '\"' << ": ";
                write(it->second, stream, indent, pos + indent, false);
                if(++it != o.end())
                    stream << ',';
                else
                    break;
            }
        stream << '\n' << std::string(pos, ' ') << '}';
    }
    
    inline void write(value const& v, std::ostream& stream, int indent = 4, int pos = 0, bool in_array = false) {
        if(v.is<null_t>())    { stream << "null";                                 return; };
        if(v.is<boolean_t>()) { stream << (v.boolean() ? "true" : "false");       return; };
        if(v.is<int64_t>())   { stream << v.int64();                              return; };
        if(v.is<real64_t>())  { stream << v.real64();                             return; };
        if(v.is<string_t>())  { stream << '\"' << v.string() << '\"';             return; };
        if(v.is<array_t>())   { write(v.array(), stream, indent, pos, in_array);  return; };
        if(v.is<object_t>())  { write(v.object(), stream, indent, pos, in_array); return; };
        if(v.is<empty_t>())     throw std::runtime_error("jsx::write: empty");
        write(v.to_json(), stream, indent, pos, in_array);
    }
    
    inline void write(value const& v, std::string const name) {
        std::ofstream file(name.c_str());  write(v, file);  file.close();
    }
    
}


#endif
