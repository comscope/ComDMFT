#ifndef JSX
#define JSX

#include <iostream>
#include <cstring>
#include <vector>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cerrno>
#include <utility>
#include <limits>
#include <map>

namespace jsx {
    
    struct TypeId {        // Saw this way of dynamic type checking somewhere on internet ... but can't find where anymore ... very nice solution !!
        using FuncPtr = TypeId(*)();
        
        TypeId() = delete;
        TypeId(FuncPtr ptr) : ptr_(ptr) {};
    private:
        FuncPtr ptr_;
        
        friend bool operator!=(TypeId const&, TypeId const&);
    };
    
    inline bool operator!=(TypeId const& lhs, TypeId const& rhs) {
        return lhs.ptr_ != rhs.ptr_;
    }
    
    template<typename T> TypeId get_type_id() { return &get_type_id<T>;};
    
    struct value;
    
    struct x_base {
        virtual void write(value&) const = 0;
        virtual x_base* copy() const = 0;
        virtual TypeId type_id() const = 0;
        virtual void* ptr() = 0;
        virtual ~x_base() {};
    };
    
    template<typename T>
    struct x_type : x_base {
        x_type() {};
        x_type(T const& t) : t_(t) {};
        x_type(T&& t) : t_(std::move(t)) {};
        void write(value& arg) const { t_.write(arg);};
        x_base* copy() const { return new x_type(t_);};
        TypeId type_id() const { return get_type_id<T>();};
        void* ptr() { return &t_;};
    private:
        T t_;
    };
    
    struct null {};
    typedef std::vector<value> array;
    typedef std::map<std::string, value> object;
    
    enum class type : uint8_t { empty, null, boolean, int64, real64, string, array, object, x };
    
    char const * const type_names[] = { "empty", "null", "boolean", "int64", "real64", "string", "array", "object", "x" };
    
    struct value {
        value() : type_(::jsx::type::empty) {};
        value(value const& arg) : type_(arg.type_) { copy(arg);};
        value(value&& arg) noexcept : type_(arg.type_), data_(arg.data_) { arg.type_ = ::jsx::type::empty;};
        value& operator=(value const& arg) { clear(); type_ = arg.type_; copy(arg); return *this;};
        value& operator=(value&& arg) noexcept { clear(); type_ = arg.type_; data_ = arg.data_; arg.type_ = ::jsx::type::empty; return *this;};
        
        /*explicit*/ value(::jsx::null const& arg) : type_(::jsx::type::null) {};
        /*explicit*/ value(bool const& arg) : type_(::jsx::type::boolean) { data_.boolean = arg;};
        /*explicit*/ value(std::int64_t const& arg) : type_(::jsx::type::int64) { data_.int64   = arg;};
        /*explicit*/ value(double const& arg) : type_(::jsx::type::real64) { data_.real64  = arg;};
        /*explicit*/ value(std::string const& arg) : type_(::jsx::type::string) { data_.string  = new std::string(arg);};
        /*explicit*/ value(::jsx::array const& arg) : type_(::jsx::type::array) { data_.array   = new ::jsx::array(arg);};
        /*explicit*/ value(::jsx::object const& arg) : type_(::jsx::type::object) { data_.object  = new ::jsx::object(arg);};
        
        /*explicit*/ value(std::string&& arg) : type_(::jsx::type::string) { data_.string  = new std::string(std::move(arg));};
        /*explicit*/ value(::jsx::array&& arg) : type_(::jsx::type::array) { data_.array   = new ::jsx::array(std::move(arg));};
        /*explicit*/ value(::jsx::object&& arg) : type_(::jsx::type::object) { data_.object  = new ::jsx::object(std::move(arg));};
        
        template <typename T, void(std::decay<T>::type::*)(value&) const = &std::decay<T>::type::write>
        /*explicit*/ value(T&& arg) : type_(::jsx::type::x) { data_.x = new x_type<typename std::decay<T>::type>(std::forward<T>(arg));}
        
        value& operator=(::jsx::null const& arg) { clear(); type_ = ::jsx::type::null; return *this;};
        value& operator=(bool const& arg) { clear(); type_ = ::jsx::type::boolean; data_.boolean = arg; return *this;};
        value& operator=(std::int64_t const& arg) { clear(); type_ = ::jsx::type::int64; data_.int64 = arg; return *this;};
        value& operator=(double const& arg) { clear(); type_ = ::jsx::type::real64; data_.real64 = arg; return *this;};
        value& operator=(std::string const& arg) { clear(); type_ = ::jsx::type::string; data_.string = new std::string(arg); return *this;};
        value& operator=(::jsx::array const& arg) { clear(); type_ = ::jsx::type::array; data_.array = new ::jsx::array(arg); return *this;};
        value& operator=(::jsx::object const& arg) { clear(); type_ = ::jsx::type::object; data_.object = new ::jsx::object(arg); return *this;};
        
        value& operator=(std::string&& arg) { clear(); type_ = ::jsx::type::string; data_.string = new std::string(std::move(arg)); return *this;};
        value& operator=(::jsx::array&& arg) { clear(); type_ = ::jsx::type::array; data_.array = new ::jsx::array(std::move(arg)); return *this;};
        value& operator=(::jsx::object&& arg) { clear(); type_ = ::jsx::type::object; data_.object = new ::jsx::object(std::move(arg)); return *this;};
        
        template <typename T, void(std::decay<T>::type::*)(value&) const = &std::decay<T>::type::write>
        value& operator=(T&& arg) { clear(); type_ = ::jsx::type::x; data_.x = new x_type<typename std::decay<T>::type>(std::forward<T>(arg)); return *this;}

        //provisorisch ?
        template<typename InputIt>
        value(InputIt begin, InputIt end) : type_(::jsx::type::array) {
            data_.array = new ::jsx::array(end - begin);
            InputIt source = begin; auto dest = array().begin();
            while(source != end) *dest++ = *source++;
        }
        
        ::jsx::type type() const { return type_;};
        
        bool boolean() const { assert_type(::jsx::type::boolean); return data_.boolean;};
        std::int64_t int64() const {
            if(type_ == ::jsx::type::real64 && static_cast<double>(static_cast<std::int64_t>(data_.real64)) == data_.real64) return data_.real64;
            assert_type(::jsx::type::int64); return data_.int64;
        };
        double real64() const {
            if(type_ == ::jsx::type::int64 && static_cast<std::int64_t>(static_cast<double>(data_.int64)) == data_.int64) return data_.int64;
            assert_type(::jsx::type::real64); return data_.real64;
        };
        std::string& string() { assert_type(::jsx::type::string); return *data_.string;};
        std::string const& string() const { assert_type(::jsx::type::string); return *data_.string;};
        ::jsx::array& array() { assert_type(::jsx::type::array); return *data_.array;};
        ::jsx::array const& array() const { assert_type(::jsx::type::array);  return *data_.array;};
        ::jsx::object& object() { assert_type(::jsx::type::object); return *data_.object;};
        ::jsx::object const& object() const { assert_type(::jsx::type::object); return *data_.object;};
        
        template <typename T, void(T::*)(value&) const = &T::write>   //google SFINAE if this looks strange to you
        T& at() {
            assert_type(::jsx::type::x);
            if(data_.x->type_id() != get_type_id<T>()) throw std::runtime_error("jsx: wrong x type");
            return *static_cast<T*>(data_.x->ptr());
        }
        
        template <typename T, void(T::*)(value&) const = &T::write>
        T const& at() const {
            assert_type(::jsx::type::x);
            if(data_.x->type_id() != get_type_id<T>()) throw std::runtime_error("jsx: wrong x type");
            return *static_cast<T const*>(data_.x->ptr());
        }
        
        void reset() { clear(); type_ = ::jsx::type::empty;}
        
        //eleganz vo arsch vo chamel aber praktisch
        std::size_t size() const {
            if(type_ == ::jsx::type::array) return array().size();
            if(type_ == ::jsx::type::object) return object().size();

            throw std::runtime_error("jsx: invalid type for size request.");
        }
        value& operator[](std::size_t index) {
            return array()[index];
        };
        value& operator()(std::size_t index) {
            if(!(index < array().size())) throw std::runtime_error("jsx: invalid array index.");
            return array().at(index);
        };
        value const& operator()(std::size_t index) const {
            if(!(index < array().size())) throw std::runtime_error("jsx: invalid array index.");
            return array().at(index);
        };
        
        std::size_t is(std::string const& key) const {
            assert_type(::jsx::type::object); return object().count(key);
        };
        value& operator[](std::string const& key) {
            if(type_ != ::jsx::type::object) *this = ::jsx::object();
            return object()[key];
        }
        value& operator()(std::string const& key) {
            if(!is(key)) throw std::runtime_error("jsx: key \"" + key + "\" not found.");
            return object().at(key);
        };
        value const& operator()(std::string const& key) const {
            if(!is(key)) throw std::runtime_error("jsx: key \"" + key + "\" not found.");
            return object().at(key);
        };
        
        void to_json(value& arg) const {
            if(type_ == ::jsx::type::x)
                data_.x->write(arg);
            else
                throw std::runtime_error("jsx::value::write: not x type.");
        };
        
        ~value() {
            clear();
        };
    private:
        ::jsx::type type_;
        
        union variant {
            bool boolean;
            std::int64_t int64;
            double real64;
            std::string* string;
            ::jsx::array* array;
            ::jsx::object* object;
            x_base* x;
        } data_;

        void assert_type(::jsx::type type) const {
            if(type != type_) throw std::runtime_error("jsx: found " +
                                                       std::string(type_names[static_cast<std::uint8_t>(type_)]) +
                                                       " type instead of requested " +
                                                       std::string(type_names[static_cast<std::uint8_t>(type)]) +
                                                       " type");
        }
        
        void clear() {
            switch(type_) {
                case ::jsx::type::string: delete data_.string; break;
                case ::jsx::type::array: delete data_.array; break;
                case ::jsx::type::object: delete data_.object; break;
                case ::jsx::type::x: delete data_.x; break;
                default: break;
            }
        }
        
        void copy(value const& arg) {
            if(type_ < ::jsx::type::string)
                data_ = arg.data_;
            else switch(arg.type_) {
                case ::jsx::type::string: data_.string = new std::string(*arg.data_.string); break;
                case ::jsx::type::array: data_.array = new ::jsx::array(*arg.data_.array); break;
                case ::jsx::type::object: data_.object = new ::jsx::object(*arg.data_.object); break;
                case ::jsx::type::x: data_.x = arg.data_.x->copy();
                default: break;
            }
        }
    };
    
    template<typename T> T& at(value& v) {
        if(v.type() != type::x) {
            T t; t.read(v); v = std::move(t);
        }
        return v.at<T>();
    };
    template<typename T> T const& at(value const& v) {
        return v.at<T>();
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
        if(*end != '\"') throw std::runtime_error("jsx parser: error while reading string");
        s.clear(); s.append(it, end);
        
        return ++end;
    }
    
    inline char const* parse(char const* it, array& a) {
        it = skip_white_space(++it);
        if(*it == ']') return ++it;
        
        while(1) {
            a.push_back(value());
            it = parse(it, a.back());
            if(*it != ',') break;
            ++it;
        };
        
        if(*it != ']') throw std::runtime_error("jsx parser: error while reading array");
        return ++it;
    }
    
    inline char const* parse(char const* it, object& o) {
        it = skip_white_space(++it);
        if(*it == '}') return ++it;
        
        std::string key;
        while(1) {
            if(*it != '\"') throw std::runtime_error("jsx parser: object key not valid");
            it = parse(it, key);
            
            it = skip_white_space(it);
            if(*it != ':') throw std::runtime_error("jsx parser: no object value found");
            it = parse(++it, o[key]);
            
            if(*it != ',') break;
            it = skip_white_space(++it);
        };
        
        if(*it != '}') throw std::runtime_error("jsx parser: error while reading object");
        return ++it;
    }
    
    const long long int lli_int64_max = std::numeric_limits<std::int64_t>::max();
    const long long int lli_int64_min = std::numeric_limits<std::int64_t>::min();
    
    inline char const* parse(char const* it, value& v) {
        it = skip_white_space(it);
        
        switch(*it) {
            case '{':
                v = object();
                it = parse(it, v.object());
                break;
            case '[':
                v = array();
                it = parse(it, v.array());
                break;
            case '\"':
                v = std::string();
                it = parse(it, v.string());
                break;
            case 'n':
                if(*++it == 'u' && *++it == 'l' && *++it == 'l') v = null(); else throw std::runtime_error("jsx parser: error while reading null");
                ++it; break;
            case 't':
                if(*++it == 'r' && *++it == 'u' && *++it == 'e') v = true; else throw std::runtime_error("jsx parser: error while reading true");
                ++it; break;
            case 'f':
                if(*++it == 'a' && *++it == 'l' && *++it == 's' && *++it == 'e') v = false; else throw std::runtime_error("jsx parser: error while reading true");
                ++it; break;
            default:
                errno = 0; char* lli_end; long long int lli = strtoll(it, &lli_end, 10);
                int int64_fail = errno != 0 || lli_end == it || lli_int64_max < lli || lli < lli_int64_min;
                
                errno = 0; char* d_end; double d = strtod(it, &d_end);
                int d_fail = errno != 0 || d_end == it;
                
                if(!int64_fail && lli_end == d_end) {
                    v = static_cast<std::int64_t>(lli); it = lli_end; break;
                }
                
                if(!d_fail) {
                    v = d; it = d_end; break;
                }
                //throw because empty ?
        }
        
        return skip_white_space(it);
    }
    
    inline void read(std::istream& stream, value& v) {
        std::string buffer;
        
        stream.seekg(0, std::ios::end);
        int size = stream.tellg();
        stream.seekg(0, std::ios::beg);
        
        buffer.resize(size + 1);
        stream.read(&buffer.front(), size);
        buffer[size] = '\0';

        if(*parse(&buffer.front(), v) != '\0') throw std::runtime_error("jsx parser: stream contains more ...");
    }
    
    inline void write(value const&, std::ostream&, int, int, bool);

    inline void write(array const& a, std::ostream& stream, int indent = 4, int pos = 0, bool in_array = false) {
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
    
    inline void write(object const& o, std::ostream& stream, int indent = 4, int pos = 0, bool in_array = false) {
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
        switch(v.type()) {
            case type::null: stream << "null"; break;
            case type::boolean: stream << (v.boolean() ? "true" : "false"); break;
            case type::int64: stream << v.int64(); break;
            case type::real64: stream << v.real64(); break;
            case type::string: stream << '\"' << v.string() << '\"'; break;
            case type::array: write(v.array(), stream, indent, pos, in_array); break;
            case type::object: write(v.object(), stream, indent, pos, in_array); break;
            case type::x: { value temp; v.to_json(temp); write(temp, stream, indent, pos, in_array); break; }
            case type::empty: throw std::runtime_error("jsx::write: empty"); break;
            default: throw std::runtime_error("jsx::write: fatal error"); break;
        }
    }
}

#endif //JSX
