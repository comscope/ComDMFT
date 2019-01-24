      subroutine pol_0_spin_spin(isp1,isp2,p)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: isp1,isp2
      real*8, intent(inout) :: p(n_pbtot,n_pbtot,ndim3_tau,
     &                           maxval(ndim3_k))
	integer :: ind_tau,ierr,idum
	real*8, allocatable :: p_r1(:,:,:),p_r2(:,:,:)
	do ind_tau=1,ndim3_tau
c	  
c	---------- MT_MT Part ------------------------------------
c	  	    
	  allocate(p_r1(n_pbmtm,n_pbmtm,nrr_max))
	  call p_mm_r(p_r1,ind_tau,isp1,isp2)
        if(ind_tau==1) 
     &      call timel('** P_R_MM finished *****************')
c ---------------------------------------------------------------------
	  allocate(p_r2(n_pbmt,n_pbmt,ndim3_kk(me3_k+1)),stat=ierr)
        if(ierr/=0) 
     &        call error_alloc('               POL_0_SPIN_SPIN',
     &                         '                          P_R2',
     &         3,n_pbmt,n_pbmt,ndim3_kk(me3_k+1),idum,idum,idum,idum,8)
	  call sym_boson_rs(p_r1,p_r2)
	  deallocate(p_r1)
	  call p_mm_q_from_r(p_r2,ind_tau,p)
	  deallocate(p_r2)
        if(ind_tau==1) 
     &      call timel('** P_MM_Q_FROM_R finished **********')
c	  
c	---------- INT_MT Part ------------------------------------
c	  
	  allocate(p_r2(nr_red,n_pbmt,ndim3_kk(me3_k+1)),stat=ierr)
        if(ierr/=0) 
     &        call error_alloc('               POL_0_SPIN_SPIN',
     &                         '                          P_R2',
     &         3,nr_red,n_pbmt,ndim3_kk(me3_k+1),idum,idum,idum,idum,8)
	  p_r2=0.d0
	  call p_im_r(ind_tau,isp1,isp2,p_r2)
        if(ind_tau==1) 
     &      call timel('** P_R_IM finished *****************')
c ---------------------------------------------------------------------
	  call p_im_q_from_r(p_r2,ind_tau,p)
	  deallocate(p_r2)
        if(ind_tau==1) 
     &      call timel('** P_IM_Q_FROM_R finished **********')
c	  
c	---------- INT_INT Part ------------------------------------
c	  
        call p_ii_solid(ind_tau,isp1,isp2,p)
        if(ind_tau==1)
     &    call timel('** P_II_SOLID finished *************')
	enddo   !! over ind_tau
      end
