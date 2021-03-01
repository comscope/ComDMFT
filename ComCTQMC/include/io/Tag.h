#ifndef INCLUDE_IO_TAG_H
#define INCLUDE_IO_TAG_H

#include <vector>
#include <complex>

#include "Vector.h"
#include "../JsonX.h"

//scheiss b64 member variable !! Kack loesig im moment

namespace io {
    
    inline void to_tagged_json(jsx::value& jArg) {
        if(jArg.is<jsx::object_t>()) {
            for(auto& jEntry : jArg.object()) to_tagged_json(jEntry.second);
        } else if(jArg.is<jsx::array_t>()) {
            for(auto& jEntry : jArg.array()) to_tagged_json(jEntry);
        } else if(!jArg.is_json()) {
            if(!(jArg.is<rvec>() || jArg.is<cvec>()))
                throw std::runtime_error("io::to_tagged_json: " + jArg.name() + " not allowed");
            jArg = jsx::object_t{{jArg.name(), jArg.to_json()}};
        }   
    }
    
    template<typename T>
    inline void read_tagged(jsx::value& jArg) {
        if(jArg.is(jsx::trait<T>::name())) {
            if(jArg.size() != 1)
                throw std::runtime_error("io::read_tagged: invalid format");
            jsx::value jTemp = jsx::at<T>(jArg(jsx::trait<T>::name())); jArg = std::move(jTemp);
        }
    }

    template<typename... Ts>
    inline void from_tagged_json(jsx::value& jArg) {
        if(jArg.is<jsx::object_t>()) read_tagged<rvec>(jArg);
        if(jArg.is<jsx::object_t>()) read_tagged<cvec>(jArg);
        if(jArg.is<jsx::object_t>()) {
            for(auto& jEntry : jArg.object()) from_tagged_json(jEntry.second);
        } else if(jArg.is<jsx::array_t>()) {
            for(auto& jEntry : jArg.array()) from_tagged_json(jEntry);
        }
    }

};

#endif
