      subroutine start_par_solid
      use comwann_mod
      implicit none
      integer :: ii
      allocate(ndim_k(nproc))
      allocate(n_mpi_k(nproc))
      allocate(ndim_kk(nproc))
      allocate(n_mpi_kk(nproc))
c     ---- Memory distribution for K-points ---------------------------------
      call size_shift_par(npnt,nproc,ndim_k,n_mpi_k)
c     ---- Memory distribution for K_All-points -----------------------------
      call size_shift_par(nqdiv,nproc,ndim_kk,n_mpi_kk)

      if (maswrk) then                                               
        write(993,*)                                                 
        write(993,'(a)') 'k point parallelization information'       
        write(993,'(a)')                                             
     $    'nproc, ndim_k, n_mpi_k, ndim_kk, n_mpi_kk'                
        do ii=1, nproc                                               
          write(993, '(5(i5, 2x))') ii, ndim_k(ii), n_mpi_k(ii),     
     $      ndim_kk(ii), n_mpi_kk(ii)                                
        enddo                                                        
        
        


      endif               
      
      end
