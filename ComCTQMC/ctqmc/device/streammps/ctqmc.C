#include "Algebra.h"
#include "../../host/Algebra.h"

#include "../../include/MonteCarlo.h"

ut::Beta ut::beta;

int main(int argc, char** argv)
{
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#endif
    try {
        if(argc != 2) throw std::runtime_error("IS: Wrong number of input parameters!");
        
        std::time_t time;
        
        mpi::cout = mpi::cout_mode::one;
        
        mpi::cout << "Start task at " << std::ctime(&(time = std::time(nullptr))) << std::endl;
        
        jsx::value jParams = mpi::read(std::string(argv[1]) + ".json");
        std::size_t const streamsPerProcess  = 1;
        std::size_t const processesPerDevice = jParams("sim per device").int64();
        
        std::vector<char> pciIds; int const deviceCount = imp::get_pci_ids(pciIds);
        std::vector<char> nodeNames = mpi::processor_name();
        
        std::cout << "Rank " << mpi::rank() << " is sitting on node " << nodeNames.data() << " and sees " << deviceCount << " device(s):";
        for(int id = 0; id < deviceCount; ++id) std::cout << " " << &pciIds[id*imp::pci_id_size()];
        std::cout << std::endl;
        
        mpi::gather(pciIds, mpi::master); mpi::gather(nodeNames, mpi::master);
        
        std::vector<char> pciId; std::vector<std::int64_t> mcIds;
        if(mpi::rank() == mpi::master) {
            pciId.resize(mpi::number_of_workers()*imp::pci_id_size(), '\0');
            
            std::map<std::string, std::map<std::string, std::size_t>> guard; std::int64_t mcId = 0;
            for(int rank = 0; rank < mpi::number_of_workers(); ++rank) {
                std::string const nodeName(&nodeNames[rank*mpi::processor_name_size()], mpi::processor_name_size());
                
                for(unsigned int id = 0; id < deviceCount; ++id) {
                    std::string const currentPciId(&pciIds[(deviceCount*rank + id)*imp::pci_id_size()], imp::pci_id_size());
                    
                    if(!guard[nodeName].count(currentPciId)) guard[nodeName][currentPciId] = 0;
                    
                    if(guard[nodeName][currentPciId] < processesPerDevice) {
                        ++guard[nodeName][currentPciId];
                        std::copy(currentPciId.begin(), currentPciId.end(), pciId.begin() + rank*imp::pci_id_size());
                        break;
                    }
                }
                
                mcIds.push_back(mcId);
                mcId += pciId[rank*imp::pci_id_size()] != '\0' ? streamsPerProcess : 1;
                mcIds.push_back(mcId);
            }
        }
        mpi::scatter(pciId, imp::pci_id_size(), mpi::master);  mpi::scatter(mcIds, 2, mpi::master);
        
        jsx::value jSimulation = jsx::array_t(mcIds.back() - mcIds.front());
        for(auto mcId = mcIds.front(); mcId < mcIds.back(); ++mcId)
            jSimulation(mcId - mcIds.front()) = jsx::object_t{
                { "id", mcId },
                { "config", jsx::read("config_" + std::to_string(mcId) + ".json", jsx::object_t()) }
            };
        
        std::int64_t mode = pciId[0] != '\0' ? 1 : 0;
        if(mode) {
            std::cout << "Rank " << mpi::rank() << " gets simulations [" << mcIds.front() << ", " << mcIds.back() << ") and uses device " << pciId.data() <<  std::endl;
            
            imp::init_device(pciId, processesPerDevice);
            mc::montecarlo<imp::Device, double>(jParams, jSimulation);
            mc::statistics<double>(jParams, jSimulation);
            imp::release_device();
        } else {
            std::cout << "Rank " << mpi::rank() << " gets simulations [" << mcIds.front() << ", " << mcIds.back() << ") and uses host" << std::endl;
            
            mc::montecarlo<imp::Host, double>(jParams, jSimulation);
            mc::statistics<double>(jParams, jSimulation);
        }
        mpi::reduce<mpi::op::sum>(mode, mpi::master);
        jSimulation["info"]["number of GPUs"] = jsx::int64_t(mode);
        
        std::vector<std::size_t> number = { jSimulation("configs").size() };  mpi::gather(number, mpi::master);
        
        mcIds.clear();
        if(mpi::rank() == mpi::master) {
            std::size_t mcId = 0;
            for(int rank = 0; rank < mpi::number_of_workers(); ++rank) {
                mcIds.push_back(mcId);
                mcId += number[rank];
                mcIds.push_back(mcId);
            }
        }
        mpi::scatter(mcIds, 2, mpi::master);
        
        for(std::size_t mcId = mcIds.front(); mcId < mcIds.back(); ++mcId) {
            std::ofstream file("config_" + std::to_string(mcId) + ".json");
            jsx::write(jSimulation("configs")(mcId - mcIds.front()), file);
            file.close();
        }
        
        mpi::write(jSimulation("measurements"), std::string(argv[1]) + ".meas.json");
        mpi::write(jSimulation("info"),         std::string(argv[1]) + ".info.json");
        
        mpi::cout << "Task of worker finished at " << std::asctime(std::localtime(&(time = std::time(nullptr)))) << std::endl;
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
    return 0;
}












