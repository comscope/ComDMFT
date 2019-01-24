#include "TraceAlgebraHost.h"

#include "../include/MarkovChain.h"
#include "../include/MonteCarlo.h"

int main(int argc, char** argv)
{
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#endif
    try {
        if(argc != 2) throw std::runtime_error("IS: Wrong number of input parameters!");
        
        std::time_t time;
        
        std::cout << "Start task at " << std::asctime(std::localtime(&(time = std::time(nullptr)))) << std::endl << std::endl;
        
        typedef tr::Comm Dummy;
        typedef ma::MarkovChain<gr::BinMoments, hy::Simple> MarkovChain;
        
        jsx::value jSimulation; mpi::read(std::string(argv[1]) + ".json", jSimulation["Parameters"]);
   
        std::vector<MarkovChain*>  markovChains(1);
        markovChains[0] = new MarkovChain(jSimulation("Parameters"), mpi::rank());

        mc::MonteCarlo<Dummy>(markovChains, jSimulation);

        delete markovChains[0];

        mpi::write(jSimulation, std::string(argv[1]) + ".meas.json");

        std::cout << "Task of worker finished at " << std::asctime(std::localtime(&(time = std::time(nullptr)))) << std::endl;
    }
    catch (ut::out_of_memory error) {
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
    
    //std::cout << tr::counter << std::endl;
    
    return 0;
}












