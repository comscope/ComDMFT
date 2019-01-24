cmake_minimum_required(VERSION 3.10.3)

# Define the project
project(ComLowH VERSION 0.0.0 LANGUAGES Fortran)

set(ComLowH_DEPENDENCIES mpi hdf5)
set(ComLowH_INCLUDES comlowh_mod.mod)

set(ComLowH_SRCS comlowh_mod.f
                 comlowh.f)

add_executable(ComLowH ${ComLowH_SRCS})

include_directories(${MPI_INCLUDE_DIR})
target_link_libraries(ComLowH ${MPI_Fortran_LIBRARIES})

include_directories(${HDF5_INCLUDE_DIR})

target_link_libraries(ComLowH ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

if(HDF5_FOUND)
  target_link_libraries(ComLowH ${HDF5_Fortran_LIBRARY_hdf5_fortran}
                                ${HDF5_Fortran_LIBRARY_hdf5}
                                ${HDF5_Fortran_LIBRARY_dl})

  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DMPI -DUSE_HDF5")
else()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DMPI")
endif()

install(TARGETS ComLowH RUNTIME DESTINATION bin)
