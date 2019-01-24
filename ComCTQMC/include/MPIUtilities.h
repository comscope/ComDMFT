#ifndef MPIUTILITIES
#define MPIUTILITIES

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstdint>

#include "JsonX.h"

namespace mpi {		
	//--------------------------------------------------schö tö tiä tü mö tiä par la barbischätöööötötötötö-------------------------------------------------------
    std::int64_t const master = 0;
	
    inline int rank() {
		int temp = 0;
		
#ifdef HAVE_MPI		
		MPI_Comm_rank(MPI_COMM_WORLD, &temp);
#endif
		
		return temp;
	}
	
    inline int number_of_workers() {
		int temp = 1;
		
#ifdef HAVE_MPI
		MPI_Comm_size(MPI_COMM_WORLD, &temp);
#endif
		
		return temp;
	}
	
	inline void read(std::string name, jsx::value& value) {
		int size; std::string buffer;
		
		if(rank() == master) {
			std::ifstream file(name.c_str()); 
			if(!file) throw std::runtime_error("mpi: file " + name + " not found !");
			
			file.seekg(0, std::ios::end);
			size = file.tellg();
			file.seekg(0, std::ios::beg);
			
            buffer.resize(size + 1);
            file.read(&buffer.front(), size);
            buffer[size] = '\0'; ++size;
            
            file.close();
        }
		
#ifdef HAVE_MPI
		MPI_Bcast(&size, 1, MPI_INT, master, MPI_COMM_WORLD);
		
		if(rank() != master) buffer.resize(size);
		
		MPI_Bcast(&buffer.front(), size, MPI_CHAR, master, MPI_COMM_WORLD);
#endif
		
        if(*jsx::parse(&buffer.front(), value) != '\0') throw std::runtime_error("json: file contains more ...");
	}	
	
    template<typename T>
	inline void write(T const& t, std::string name) {
		if(rank() == master) {
			std::ofstream file(name.c_str());
			jsx::write(t, file);
			file.close();
        } else {
            std::ostringstream dummy;
            jsx::write(t, dummy);
        }
        
#ifdef HAVE_MPI
		MPI_Barrier(MPI_COMM_WORLD);    //Why the hell this ?!?!??!
#endif
	}
    

    inline void broad_cast(jsx::value& value) {
#ifdef HAVE_MPI
        int size; std::string buffer;
        
        if(rank() == master) {
            std::ostringstream stream;
            jsx::write(value, stream);
            buffer = stream.str();
            buffer.push_back('\0');
            size = buffer.size();
        }
        
        MPI_Bcast(&size, 1, MPI_INT, master, MPI_COMM_WORLD);
        
        if(rank() != master) buffer.resize(size);
        
        MPI_Bcast(&buffer.front(), size, MPI_CHAR, master, MPI_COMM_WORLD);
        
        if(*jsx::parse(&buffer.front(), value) != '\0') throw std::runtime_error("json: file contains more ...");
#endif    
    };
}


#endif
