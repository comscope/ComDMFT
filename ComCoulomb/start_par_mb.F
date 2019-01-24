      subroutine start_par_mb
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      include 'mpif.h'
      integer :: ierr
c     ---- Memory distribution for Double Parallelization  -----------------------
      if(mod(nproc,nproc_tau*nproc_k)/=0) then
        if(me==0) write(iun,*)' NPROC does not match NPROC_TAU*NPROC_K'
c        write(179,*) "calling ending in start_par_mb.f"
        call ending
      endif
      if(mod(n_tau/2+1,nproc_tau)/=0) then
        if(me==0) write(iun,*)' N_TAU does not match NPROC_TAU'
c        write(179,*) "calling ending in start_par_mb.f"
        call ending
      endif
      if(mod(n_omega+1,nproc_tau)/=0) then
        if(me==0) write(iun,*)' N_OMEGA does not match NPROC_TAU'
c        write(179,*) "calling ending in start_par_mb.f"
        call ending
      endif
      if(mod(n_nu+1,nproc_tau)/=0) then
        if(me==0) write(iun,*)' N_NU does not match NPROC_TAU'
c        write(179,*) "calling ending in start_par_mb.f"
        call ending
      endif
      ndim3_tau=(n_tau/2+1)/nproc_tau
      me3_tau=mod(me,nproc_tau)
      ndim3_omega=(n_omega+1)/nproc_tau
      ndim3_nu=(n_nu+1)/nproc_tau

      ndim3_rfreq=nint((maxrealfreq/(0.1d0/evolt*2.0d0))/nproc_tau)
      if (mod(ndim3_rfreq , 2) .ne. 0) then
        ndim3_rfreq=ndim3_rfreq+1
      endif
      n_rfreq=ndim3_rfreq*nproc_tau
      
      if (maswrk) write(993, *) 'ndim3_tau=', ndim3_tau
      if (maswrk) write(993, *) 'ndim3_omega=', ndim3_omega
      if (maswrk) write(993, *) 'ndim3_nu=', ndim3_nu
      if (maswrk) write(993, *) 'ndim3_rfreq=', ndim3_rfreq
      allocate(ndim_tau(nproc))
      allocate(ndim_omega(nproc))
      allocate(ndim_nu(nproc))
      allocate(ndim_rfreq(nproc))
      allocate(ndim_istar(nproc))
      allocate(n_mpi_tau(nproc))
      allocate(n_mpi_omega(nproc))
      allocate(n_mpi_nu(nproc))
      allocate(n_mpi_rfreq(nproc))
      allocate(n_mpi_istar(nproc))
!      call mpi_bcast(chem_pot,1,mpi_double_precision,
!     $    0,MPI_COMM_WORLD,mpierr)
!      write(179,*)'bk-start_par 1'
c     ---- Memory distribution for Omega-mesh ---------------------------------
      call size_shift_par(n_omega+1,nproc,ndim_omega,n_mpi_omega)
c     ---- Memory distribution for Tau-mesh ---------------------------------
      call size_shift_par(n_tau/2+1,nproc,ndim_tau,n_mpi_tau)
c     ---- Memory distribution for Nu-mesh ---------------------------------
      call size_shift_par(n_nu+1,nproc,ndim_nu,n_mpi_nu)
c     ---- Memory distribution for rfreq-mesh ---------------------------------
      call size_shift_par(n_rfreq,nproc,ndim_rfreq,n_mpi_rfreq)
c     ------------ Communicators ----------------------
      me3_k=mod(me/nproc_tau,nproc_k)
!      call mpi_bcast(chem_pot,1,mpi_double_precision,
!     $    0,MPI_COMM_WORLD,mpierr)
!      write(179,*)'bk-start_par 2'

      call MPI_COMM_SPLIT(MPI_COMM_WORLD,me3_tau,me,comm_tau,ierr)
c$$$      call comm_rank(comm_tau,me3_k)
      call mpi_comm_rank(comm_tau,me3_k, mpierr)
      call MPI_COMM_SPLIT(MPI_COMM_WORLD,me3_k,me,comm_k,ierr)
c$$$      call comm_rank(comm_k,me3_tau)
      call mpi_comm_rank(comm_k,me3_tau, mpierr)  
      
!      call mpi_bcast(chem_pot,1,mpi_double_precision,
!     $    0,MPI_COMM_WORLD,mpierr)
!      write(179,*)'bk-start_par 3'

c      write(179,*)'bk-start_par end'    
      end
