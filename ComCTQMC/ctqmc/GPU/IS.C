#include "TraceAlgebraHost.h"

#define SWITCH(arg) arg##Host
#include "../include/MarkovChain.h"
#undef SWITCH

#undef TRACEALGEBRA     //huiuiui schön isch andersch, aber dä ganzi mischt bis ufe zur TraceAlgebra vertemplate isch au kacke.
#undef TRACEELEMENT
#undef TRACEUTILITIES
#undef DYNAMICTRACE
#undef TRACE
#undef MARKOVCHAIN

#include "TraceAlgebraDevice.h"

#define SWITCH(arg) arg##Device
#include "../include/MarkovChain.h"
#undef SWITCH

#include "../include/MonteCarlo.h"


int main(int argc, char** argv)
{
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#endif
    try {
        if(argc != 2) throw std::runtime_error("IS: Wrong number of input parameters!");
        
        std::time_t time;
        
        std::cout << "Start task at " << std::ctime(&(time = std::time(NULL))) << std::endl;
        
        Ut::Simulation simulation(argv[1]);
        Json::Value const& jParams = simulation.params();
        
        int const core_per_device = jParams("CORE_PER_DEVICE").int64();
        int const device_per_node = jParams("DEVICE_PER_NODE").int64();
        int const stream_per_device = jParams("STREAM_PER_DEVICE").int64();
        int const stream_minus_rank = ((mpi::rank() + core_per_device - 1)/core_per_device)*(stream_per_device - 1);
        
        if(mpi::rank() % core_per_device == 0) {
            typedef MaDevice::Comm Comm; typedef MaDevice::MarkovChain<Green::BinMoments, Hyb::Simple> MarkovChain;

            Comm::init((mpi::rank() / core_per_device) % device_per_node, stream_per_device, jParams("MEMORY_PER_DEVICE").real64(), 4096);
            
            MarkovChain::init(jParams);
        
            std::vector<MarkovChain*>  markovChains(stream_per_device);
            for(std::size_t s = 0; s < markovChains.size(); ++s)
                markovChains[s] = new MarkovChain(jParams, mpi::rank() + s + stream_minus_rank);
            
            MC::MonteCarlo<Comm>(markovChains, simulation);
            
            for(std::size_t s = 0; s < markovChains.size(); ++s)
                delete markovChains[s];
           
            MarkovChain::release();
            Comm::release();
        } else {
            typedef MaHost::Comm Dummy; typedef MaHost::MarkovChain<Green::BinMoments, Hyb::Simple> MarkovChain;
            
            MarkovChain::init(jParams);
            
            std::vector<MarkovChain*>  markovChains(1);
            markovChains[0] = new MarkovChain(jParams, mpi::rank() + stream_minus_rank);
            
            MC::MonteCarlo<Dummy>(markovChains, simulation);
            
            delete markovChains[0];
            
            MarkovChain::release();
        }
        
        std::cout << "Task of worker finished at " << std::ctime(&(time = std::time(NULL))) << std::endl;
    }
    catch (Ut::out_of_memory error) {
        std::cerr << "Out of memory" << std::endl;
        
#ifdef HAVE_MPI
        MPI_Abort(MPI_COMM_WORLD, -1);
#endif
        return -1;
    }
    catch (std::exception& exc) {
        std::cerr << exc.what() << " ( Thrown from worker " << mpi::rank() << " )" << std::endl;
        
#ifdef HAVE_MPI
        MPI_Abort(MPI_COMM_WORLD, -1);
#endif
        return -1;
    }
    catch (...) {
        std::cerr << "Fatal Error: Unknown Exception! ( Thrown from worker " << mpi::rank() << " )" << std::endl;
        
#ifdef HAVE_MPI
        MPI_Abort(MPI_COMM_WORLD, -2);
#endif
        return -2;
    }
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}












