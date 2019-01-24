      Program QFT
	use manager_mod
	use parallel_mod
      use files_mod
#ifdef USE_HDF5
      use hdf5io_mod
#endif
c***********************************************************
!              The main module                             *
!***********************************************************
      implicit none

      call beging
      call debugger
	mem_using=0.d0
	mem_max=0.d0
	mem_where='qft'
#ifdef USE_HDF5
      call finit()
#endif
      call init
      if(maswrk) write(iun,*)' Number of processors is ',nproc
	call start_par_mb
        call tau_frequency_meshes


! --------------------------------------------------------------
      call solid
!     ***********************************************************
#ifdef USE_HDF5
c      call fclose(fid_rst)
c      call fclose(fid_dat)
      call ffinal()
#endif
      call timel("**** PROGRAM FINISHED **************")
      call pend
	end
