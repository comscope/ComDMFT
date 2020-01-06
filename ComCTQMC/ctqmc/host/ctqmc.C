#include "AlgebraHost.h"

#include "../include/MonteCarlo.h"

ut::Beta ut::beta;

int main(int argc, char** argv)
{
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#endif
    try {
        if(argc != 2) throw std::runtime_error("ctqmc: Wrong number of input parameters!");
        
        std::time_t time;
        
        mpi::cout = mpi::cout_mode::one;
        
        mpi::cout << "Start task at " << std::asctime(std::localtime(&(time = std::time(nullptr)))) << std::endl << std::endl;

        jsx::value jParams = mpi::read(std::string(argv[1]) + ".json");
        
        jsx::value jSimulation = jsx::array_t(1);
        jSimulation(0) = jsx::int64_t(mpi::rank());
        
        if(jParams.is("complex hybridisation") ? jParams("complex hybridisation").boolean() : false)
            mc::MonteCarlo<imp::AllocHost, ut::complex>(jParams, jSimulation);
        else
            mc::MonteCarlo<imp::AllocHost, double>(jParams, jSimulation);
        
        mpi::write(jSimulation("measurements"), std::string(argv[1]) + ".meas.json");
        mpi::write(jSimulation("info"),         std::string(argv[1]) + ".info.json");
        
        std::ofstream file(("config_" + std::to_string(mpi::rank()) + ".json").c_str());
        jsx::write(jSimulation("configs")(0), file);
        file.close();

        mpi::cout << "Task of worker finished at " << std::asctime(std::localtime(&(time = std::time(nullptr)))) << std::endl;
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












