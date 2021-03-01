#ifndef CTQMC_INCLUDE_UPDATES_INCLUDE_SURVIVING_H
#define CTQMC_INCLUDE_UPDATES_INCLUDE_SURVIVING_H

#include "../../config/Worms.h"

namespace upd {
    
    //---------------------------------------------------------------  forward declaration  ----------------------------------------------------------------------------------
    
    template<typename Value> void surviving_map(imp::itf::Operators<Value> const&, std::vector<int>&, cfg::Flavor const&);
    template<typename Value> void surviving_map(imp::itf::Operators<Value> const&, std::vector<int>&, std::tuple<cfg::Flavor, cfg::Flavor> const&);
    template<typename Value, typename Arg1, typename Arg2, typename... Args> bool surviving_map(imp::itf::Operators<Value> const&, std::vector<int>, Arg1 const&, Arg2 const&, Args const&...);
    
    
    //-----------------------------------------------------------------  all permutations  -----------------------------------------------------------------------------------
    
    template<typename Value, typename Arg>
    bool surviving_perm(imp::itf::Operators<Value> const& ops, std::vector<int>& sectors, Arg const& arg) {
        surviving_map(ops, sectors, arg);
        for(int sec = 1; sec < sectors.size(); ++sec) if(sectors[sec] == sec) return true;
        return false;
    }
    
    template<typename Value, typename Arg1, typename Arg2>
    bool surviving_perm(imp::itf::Operators<Value> const& ops, std::vector<int>& sectors, Arg1 const& arg1, Arg2 const& arg2) {
        if(surviving_map(ops, sectors, arg1, arg2)) return true;
        if(surviving_map(ops, sectors, arg2, arg1)) return true;
        return false;
    }
    
    template<typename Value, typename Arg1, typename Arg2, typename Arg3>                                             // that is enough for the moment ....
    bool surviving_perm(imp::itf::Operators<Value> const& ops, std::vector<int>& sectors, Arg1 const& arg1, Arg2 const& arg2, Arg3 const& arg3) {
        if(surviving_map(ops, sectors, arg1, arg2, arg3)) return true;
        if(surviving_map(ops, sectors, arg2, arg1, arg3)) return true;
        if(surviving_map(ops, sectors, arg3, arg2, arg1)) return true;
        return false;
    }
    
    
    //--------------------------------------------------------------------  map sectors  ------------------------------------------------------------------------------------
    
    template<typename Value>
    void surviving_map(imp::itf::Operators<Value> const& ops, std::vector<int>& sectors, cfg::Flavor const& arg) {
        int const flavor = arg.flavor()%2 ? arg.flavor() - 1 : arg.flavor() + 1;                                      // take hermitian conjugate because we are going in the wrong
        for(auto& sec : sectors) if(sec != 0) sec = ops.at(flavor).map(sec).sector;                                   // direction (operator maps <- but we are going ->)
    };
    
    template<typename Value>
    void surviving_map(imp::itf::Operators<Value> const& ops, std::vector<int>& sectors, std::tuple<cfg::Flavor, cfg::Flavor> const& args) {
        surviving_map(ops, sectors, std::get<0>(args));
        surviving_map(ops, sectors, std::get<1>(args));
    };
    
    template<typename Value, typename Arg1, typename Arg2, typename... Args>
    bool surviving_map(imp::itf::Operators<Value> const& ops, std::vector<int> sectors, Arg1 const& arg1, Arg2 const& arg2, Args const&... args) {
        surviving_map(ops, sectors, arg1);
        return surviving_perm(ops, sectors, arg2, args...);
    };
    
    
    //----------------------------------------------------------------------  start  ------------------------------------------------------------------------------------------
    
    template<typename Value, typename Arg>
    bool surviving(data::Data<Value> const& data, Arg const& arg) {
        std::vector<int> sectors(data.eig().sectorNumber() + 1); std::iota(sectors.begin(), sectors.end(), 0);
        return surviving_perm(data.ops(), sectors, arg);
    };
    
    template<typename Value, typename Arg1, typename Arg2, typename... Args>
    bool surviving(data::Data<Value> const& data, Arg1 const& arg1, Arg2 const& arg2, Args const&... args) {
        std::vector<int> sectors(data.eig().sectorNumber() + 1); std::iota(sectors.begin(), sectors.end(), 0);
        return surviving_map(data.ops(), sectors, arg1, arg2, args...);                                               // no need to permute arg1 because trace is cyclic
    };
    
    
}


#endif
