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
    
    struct value;
    
    struct TypeId {        // Saw this way of dynamic type checking somewhere on internet ... but can't find where anymore ... very nice solution !!
        using FuncPtr = TypeId(*)();
        TypeId(FuncPtr val) : val_(val) {};
    private:
        FuncPtr val_;
        friend value;
    };

    template<typename T> TypeId get_type_id() { return &get_type_id<T>;};

    
    template<typename T> struct Trait {
        enum : bool { is_json = false };
        static std::string name(T const& t = T()) { return T::name();};
        static void write(T const& t, value& dest) { t.write(dest);};
    };

    struct empty_t {};                              inline std::string json_name(empty_t const&)   { return "jsx::empty_t";};
    struct null_t {};                               inline std::string json_name(null_t const&)    { return "jsx::null_t";};
    typedef bool boolean_t;                         inline std::string json_name(boolean_t const&) { return "jsx::boolean_t";};
    typedef std::int64_t int64_t;                   inline std::string json_name(int64_t const&)   { return "jsx::int64_t";};
    typedef double real64_t;                        inline std::string json_name(real64_t const&)  { return "jsx::real64_t";};
    typedef std::string string_t;                   inline std::string json_name(string_t const&)  { return "jsx::string_t";};
    typedef std::vector<value> array_t;             inline std::string json_name(array_t const&)   { return "jsx::array_t";};
    typedef std::map<std::string, value> object_t;  inline std::string json_name(object_t const&)  { return "jsx::object_t";};
    
    template<typename T> struct IsJsonTrait {
        enum : bool { is_json = true };
        static std::string name(T const& t = T()) { return json_name(t);};
        static void write(T const& t, value& dest) { throw std::runtime_error("jsx: writing json type is invalid"); };
    };
    
    template<> struct Trait<empty_t>   : IsJsonTrait<empty_t>   {};
    template<> struct Trait<null_t>    : IsJsonTrait<null_t>    {};
    template<> struct Trait<boolean_t> : IsJsonTrait<boolean_t> {};
    template<> struct Trait<int64_t>   : IsJsonTrait<int64_t>   {};
    template<> struct Trait<real64_t>  : IsJsonTrait<real64_t>  {};
    template<> struct Trait<string_t>  : IsJsonTrait<string_t>  {};
    template<> struct Trait<array_t>   : IsJsonTrait<array_t>   {};
    template<> struct Trait<object_t>  : IsJsonTrait<object_t>  {};
    

    struct Any {
        virtual std::string name() const = 0;
        virtual bool is_json() const = 0;
        virtual void write(value&) const = 0;
        virtual Any* clone() const = 0;
        virtual TypeId type_id() const = 0;
        virtual void* ptr() = 0;
        virtual ~Any() = default;
    };
    
    template<typename T>
    struct AnyImpl : Any {
        AnyImpl(T const& t) : t_(t) {};
        AnyImpl(T&& t) : t_(std::move(t)) {};
        std::string name() const { return Trait<T>::name(t_);};
        bool is_json() const { return Trait<T>::is_json;};
        void write(value& dest) const { Trait<T>::write(t_, dest);};
        Any* clone() const { return new AnyImpl(t_);};
        TypeId type_id() const { return get_type_id<T>();};
        void* ptr() { return &t_;};
    private:
        T t_;
    };

    struct value {
        value() : data_(new AnyImpl<empty_t>(empty_t())) {};
        value(value const& arg) : data_(arg.data_->clone()) {};
        value(value&& arg) noexcept : data_(arg.data_) { arg.data_ = new AnyImpl<empty_t>(empty_t());};
        value& operator=(value const& arg) { delete data_; data_ = arg.data_->clone(); return *this;};
        value& operator=(value&& arg) noexcept { delete data_; data_ = arg.data_; arg.data_ = new AnyImpl<empty_t>(empty_t()); return *this;};
        ~value() { delete data_;};
        
        template<typename T, typename std::enable_if<!std::is_same<typename std::decay<T>::type, value>::value, int>::type = 0>
        value(T&& arg) : data_(new AnyImpl<typename std::decay<T>::type>(std::forward<T>(arg))) {
        }
        template<typename T, typename std::enable_if<!std::is_same<typename std::decay<T>::type, value>::value, int>::type = 0>
        value& operator=(T&& arg) {
            delete data_; data_ = new AnyImpl<typename std::decay<T>::type>(std::forward<T>(arg)); return *this;
        }
        
        template<typename InputIt> //provisorisch ?
        value(InputIt begin, InputIt end) : data_(new AnyImpl<array_t>(array_t(end - begin))) {
            InputIt source = begin; auto dest = array().begin();
            while(source != end) *dest++ = *source++;
        }

        template<typename T>
        bool is() const {
            return data_->type_id().val_ == get_type_id<T>().val_;
        }
        
        template<typename T> T& at() { return at_impl<T>();}
        template<typename T> T const& at() const { return at_impl<T>();}
        
        boolean_t& boolean() { return at<boolean_t>();};
        boolean_t const& boolean() const { return at<boolean_t>();};
        string_t& string() { return at<string_t>();};
        string_t const& string() const { return at<string_t>();};
        array_t& array() { return at<array_t>();};
        array_t const& array() const { return at<array_t>();};
        object_t& object() { return at<object_t>();};
        object_t const& object() const { return at<object_t>();};
        
        int64_t int64() const {
            if(is<real64_t>() && static_cast<real64_t>(static_cast<int64_t>(unsafe_at_impl<real64_t>())) == unsafe_at_impl<real64_t>()) return unsafe_at_impl<real64_t>();
            return at<int64_t>();
        };
        real64_t real64() const {
            if(is<int64_t>() && static_cast<int64_t>(static_cast<real64_t>(unsafe_at_impl<int64_t>())) == unsafe_at_impl<int64_t>()) return unsafe_at_impl<int64_t>();
            return at<real64_t>();
        };
        
        std::size_t size() const { //eleganz vo arsch vo chamel aber praktisch
            if(is<array_t>()) return array().size();
            if(is<object_t>()) return object().size();

            throw std::runtime_error("jsx::value: found " + data_->name() + " type instead of jsx::array_t or jsx::object_ type for size request.");
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
            return data_->is_json();
        };

    private:
        Any* data_;
        
        template<typename T> T& unsafe_at_impl() const {
            return *static_cast<T*>(data_->ptr());
        }
        
        template<typename T>
        T& at_impl() const {
            if(!is<T>()) throw std::runtime_error("jsx::value: found " + data_->name() + " type instead of " + Trait<T>::name() + " type.");
            return unsafe_at_impl<T>();
        }
        
        friend value to_jsx(value const&);
    };

    
    template<typename T, typename std::enable_if< !Trait<T>::is_json, int>::type = 0>
    T& at(value& source) {
        if(source.is_json()) {
            T dest; dest.read(source); source = std::move(dest);
        }
        return source.at<T>();
    };
    template<typename T, typename std::enable_if< !Trait<T>::is_json, int>::type = 0>
    T const& at(value const& source) {
        return source.at<T>();
    };

    inline value to_jsx(value const& source){
        value dest; source.data_->write(dest); return dest;
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
        
        s.append(end, end + 1);
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
    
    template<typename T>
    inline void read(std::istream& stream, T& v) {
        std::string buffer;
        
        stream.seekg(0, std::ios::end);
        int size = stream.tellg();
        stream.seekg(0, std::ios::beg);
        
        buffer.resize(size + 1);
        stream.read(&buffer.front(), size);
        buffer[size] = '\0';

        if(*parse(&buffer.front(), v) != '\0') throw std::runtime_error("jsx::parser: stream contains more ...");
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
        if(v.is<null_t>())    { stream << "null"; return;};
        if(v.is<boolean_t>()) { stream << (v.boolean() ? "true" : "false"); return;};
        if(v.is<int64_t>())   { stream << v.int64(); return;};
        if(v.is<real64_t>())  { stream << v.real64(); return;};
        if(v.is<string_t>())  { stream << '\"' << v.string() << '\"'; return;};
        if(v.is<array_t>())   { write(v.array(), stream, indent, pos, in_array); return;};
        if(v.is<object_t>())  { write(v.object(), stream, indent, pos, in_array); return;};
        if(v.is<empty_t>())   { throw std::runtime_error("jsx::write: empty");};
        value temp = to_jsx(v); write(temp, stream, indent, pos, in_array);
    }
}

#endif //JSX
