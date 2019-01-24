#ifndef TRACEALGEBRAHOST
#define TRACEALGEBRAHOST

#include <list>
#include <set>
#include <stdexcept>
#include <vector>

namespace TrHost {
    struct Comm {
        static int iStream() { return iStream_;};
        static int is_ready(int iStream) { iStream_ = iStream; return 1;};
        static void launch() {};
    private:
        static int iStream_;
    };
    
    typedef double* data_ptr;
    
    inline double* get(data_ptr ptr) { return ptr;};
    
}


#endif
