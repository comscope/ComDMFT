cmake_minimum_required(VERSION 3.10.3)

# Define the project
project(ComWann VERSION 0.0.0 LANGUAGES Fortran)

set(ComWann_DEPENDENCIES mpi hdf5 wannier lapack blas)

set(ComWann_SRCS comwann_mod.F
                 hdf5io_mod.F
                 beging.F
                 bessr.F
                 cft.F
                 cross.F
                 det3.F
                 diagonalization_lapack.F
                 dif.F
                 dqdall.F
                 fac2.F
                 fft3.F
                 fifi_j_prepare.F
                 gbs_number.F
                 getlimj.F
                 init.F
                 integral_band_pw_band.F
                 integral_band_pw_band_mt.F
                 lget.F
                 output_u_f0f2f4f6.F
                 prd1.F
                 prd2.F
                 qft.F
                 read_restart_file.F
                 read_restart_file_hdf5.F
                 read_wan_bni.F
                 read_wan_bni_hdf5.F
                 read_wan_ini.F
                 rotate.F
                 rotmat_wan_proj.F
                 set_memory.F
                 size_shift_par.F
                 solid.F
                 sphharm.F
                 start_par_solid.F
                 sym_a_0.F
                 sym_z_0.F
                 tcgstb.F
                 unport.F
                 wannier_amn_from_win.F
                 wannier_interface.F
                 wannier_mmn.F
                 wan_pb_rs.F)

add_executable(ComWann ${ComWann_SRCS})

include_directories(${MPI_INCLUDE_DIR})
target_link_libraries(ComWann ${MPI_Fortran_LIBRARIES})

include_directories(${HDF5_INCLUDE_DIR})

target_link_libraries(ComWann $<TARGET_OBJECTS:libwannier>)

target_link_libraries(ComWann ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})


if(HDF5_FOUND)
  target_link_libraries(ComWann ${HDF5_Fortran_LIBRARY_hdf5_fortran}
                                ${HDF5_Fortran_LIBRARY_hdf5}
                                ${HDF5_Fortran_LIBRARY_dl})

  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DMPI -DUSE_HDF5")
else()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DMPI")
endif()

install(TARGETS ComWann RUNTIME DESTINATION bin)
