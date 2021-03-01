#ifndef CTQMC_INCLUDE_CONFIG_WORMS_H
#define CTQMC_INCLUDE_CONFIG_WORMS_H


#include <stdexcept>
#include <iostream>
#include <array>
#include <vector>

#include "Variant.h"
#include "Tuple.h"


namespace cfg {
    
    
    namespace partition {
        
        constexpr char name[] = "partition";
        
        using Worm = WormTuple<name>;
        
    }
    
    namespace green {
        
        constexpr char name[] = "green";
        
        using Worm = WormTuple<name, op, opDagg>;
        
    }
    
    namespace green_impr {
        
        constexpr char name[] = "green impr";
        
        using Worm = WormTuple<name, opBulla, opDagg>;
        
    }
    
    namespace green_imprsum {
        
        constexpr char name[] = "green imprsum";
        
        using Worm = WormTuple<name, opBullaSum, opDagg>;
        
    }
    
    namespace vertex {
        
        constexpr char name[] = "vertex";
        
        using Worm = WormTuple<name, op, opDagg, op, opDagg>;
        
    }
    
    namespace vertex_impr {
        
        constexpr char name[] = "vertex impr";
        
        using Worm = WormTuple<name, opBulla, opDagg, op, opDagg>;
        
    }
    
    namespace vertex_imprsum {
        
        constexpr char name[] = "vertex imprsum";
        
        using Worm = WormTuple<name, opBullaSum, opDagg, op, opDagg>;
        
    }
    
    namespace susc_ph {
        
        constexpr char name[] = "susc ph";
        
        using Worm = WormTuple<name, bilinearPH, bilinearPH>;
        
    }
    
    namespace susc_pp {
        
        constexpr char name[] = "susc pp";
        
        using Worm = WormTuple<name, bilinearHH, bilinearPP>;
        
    }
    
    namespace hedin_ph {
        
        constexpr char name[] = "hedin ph";
        
        using Worm = WormTuple<name, op, opDagg, bilinearPH>;
        
    }
    
    namespace hedin_ph_impr {
        
        constexpr char name[] = "hedin ph impr";
        
        using Worm = WormTuple<name, opBulla, opDagg, bilinearPH>;
        
    }
    
    namespace hedin_ph_imprsum {
        
        constexpr char name[] = "hedin ph imprsum";
        
        using Worm = WormTuple<name, opBullaSum, opDagg, bilinearPH>;
        
    }
    
    namespace hedin_pp {
        
        constexpr char name[] = "hedin pp";
        
        using Worm = WormTuple<name, op, op, bilinearPP>;
        
    }
    
    namespace hedin_pp_impr {
        
        constexpr char name[] = "hedin pp impr";
        
        using Worm = WormTuple<name, opBulla, op, bilinearPP>;
        
    }
    
    namespace hedin_pp_imprsum {
        
        constexpr char name[] = "hedin pp imprsum";
        
        using Worm = WormTuple<name, opBullaSum, op, bilinearPP>;
        
    }
    
    
    using Worm = WormVariant<
    
    partition::Worm,
    
    green::Worm,
    green_impr::Worm,
    green_imprsum::Worm,
    
    vertex::Worm,
    vertex_impr::Worm,
    vertex_imprsum::Worm,
    
    susc_ph::Worm,
    
    susc_pp::Worm,
    
    hedin_ph::Worm,
    hedin_ph_impr::Worm,
    hedin_ph_imprsum::Worm,
    
    hedin_pp::Worm,
    hedin_pp_impr::Worm,
    hedin_pp_imprsum::Worm
    
    >;
    
}

#endif
