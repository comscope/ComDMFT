      subroutine start_par_solid
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer :: ii
      allocate(ndim_k(nproc))
      allocate(n_mpi_k(nproc))
      allocate(ndim_kk(nproc))
      allocate(n_mpi_kk(nproc))
      allocate(ndim3_k(nproc_k))
      allocate(n3_mpi_k(nproc_k))
      allocate(ndim3_kk(nproc_k))
      allocate(n3_mpi_kk(nproc_k))
      allocate(ndim3_kkt(nproc_tau))
      allocate(n3_mpi_kkt(nproc_tau))
c     ---- Memory distribution for K-points ---------------------------------
      call size_shift_par(npnt,nproc,ndim_k,n_mpi_k)
      call size_shift_par(npnt,nproc_k,ndim3_k,n3_mpi_k)
c     ---- Memory distribution for K_All-points -----------------------------
      call size_shift_par(nqdiv,nproc,ndim_kk,n_mpi_kk)
      call size_shift_par(nqdiv,nproc_k,ndim3_kk,n3_mpi_kk)
      call size_shift_par(nqdiv,nproc_tau,ndim3_kkt,n3_mpi_kkt)
      allocate(ndim3_nr_red(nproc_k))
      allocate(n3_mpi_nr_red(nproc_k))
      call size_shift_par(nr_red,nproc_k,ndim3_nr_red,n3_mpi_nr_red)
      allocate(ndimt_nr_red(nproc_tau))
      allocate(nt_mpi_nr_red(nproc_tau))
      call size_shift_par(nr_red,nproc_tau,ndimt_nr_red,nt_mpi_nr_red)
c     ----------------------------------------------------------------------
      allocate(ndim_nfun(nproc_k))
      allocate(n_mpi_nfun(nproc_k))
      call size_shift_par(nfun,nproc_k,ndim_nfun,n_mpi_nfun)
      allocate(ndim_nrr(nproc_k))
      allocate(n_mpi_nrr(nproc_k))
      call size_shift_par(nrr_max,nproc_k,ndim_nrr,n_mpi_nrr)
      allocate(ndim_nrr_cor(nproc_k))
      allocate(n_mpi_nrr_cor(nproc_k))
      call size_shift_par(nrr_max_cor,nproc_k,ndim_nrr_cor,
     &  n_mpi_nrr_cor)
      allocate(ndim_omega_d(nproc_k))
      allocate(n_mpi_omega_d(nproc_k))
      call size_shift_par(n_omega_d+1,nproc_k,ndim_omega_d,
     &  n_mpi_omega_d)
      allocate(ndim_nu_d(nproc_k))
      allocate(n_mpi_nu_d(nproc_k))
      call size_shift_par(n_nu_d+1,nproc_k,ndim_nu_d,n_mpi_nu_d)



      if (maswrk) then                                               
        write(993,*)                                                 
        write(993,'(a)') 'k point parallelization information'       
        write(993,'(a)')                                             
     $    'nproc, ndim_k, n_mpi_k, ndim_kk, n_mpi_kk'                
        do ii=1, nproc                                               
          write(993, '(5(i5, 2x))') ii, ndim_k(ii), n_mpi_k(ii),     
     $      ndim_kk(ii), n_mpi_kk(ii)                                
        enddo                                                        
        
        
        write(993,'(a)')                                             
     $    'nproc_k, ndim3_k, n3_mpi_k, ndim3_kk, n3_mpi_kk'          
        do ii=1, nproc_k                                             
          write(993, '(5(i5, 2x))') ii,                              
     $      ndim3_k(ii), n3_mpi_k(ii), ndim3_kk(ii), n3_mpi_kk(ii)   
        enddo                                                        


        write(993,*)
        write(993,*)                                                 
        write(993,'(a)') 'nfun information'       
        write(993,'(a)')                                             
     $    'nproc, ndim_k, n_mpi_k, ndim_kk, n_mpi_kk'                
        do ii=1, nproc_k
          write(993, '(3(i5, 2x))') ii, ndim_nfun(ii), n_mpi_nfun(ii)
        enddo                                                        
      endif               
      
      end
