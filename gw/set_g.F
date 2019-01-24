      subroutine set_g
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      logical :: correlation
      character*8 :: way
	integer :: k,ispin,n,ind_omega,i_omega,i_tau,ind_k,i,it,ind_tau
	real*8 :: w_n,ta,cf(2)
	real*8, allocatable :: gx_tau(:,:)
	complex*16, allocatable :: gx_omega(:),tmp(:,:)
	correlation=.false.
      if(ubi==' qp'.or.ubi==' gw')
     &  correlation=.true.
	call etot_gw_0
	allocate(tmp(nbndf,nbndf))
      if(correlation) then
	  allocate(gc_omega(nbndf,nbndf,2,ndim3_omega))
	  allocate(gc_tau(nbndf,nbndf,2,ndim3_tau))
      endif
     	if(ubi==' gw') way='GW_based'
      if(ubi==' qp') way='QP_based'  
	do ispin=1,nspin
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    n=n_bnd(k,ispin)
	    allocate(gx_tau(n,2))
          call g_x_tau(ispin,k,gx_tau,betta_t,n,chem_pot)
          if(correlation) then
            if(way=='GW_based') then 
	        allocate(gx_omega(n))
	        gc_omega=0.d0
	        do ind_omega=1,ndim3_omega
	          i_omega=me3_tau*ndim3_omega+ind_omega-1
	          w_n=w_omega(i_omega)
c ------ We temporarily place SIGMA_C into TMP --------------------
		        call ferm_unpack_omega(tmp,
     &		                  sig_c_omega(1,1,1,ind_omega,ind_k,ispin),
     &                                 n,nbndf,nbndf)
		        call g_x_omega(ispin,k,gx_omega,w_n,n,chem_pot)
                call g_c_omega(0,gx_omega,tmp,n,nbndf)
			    call ferm_pack_omega(tmp,gc_omega(1,1,1,ind_omega),n,
     &			                     nbndf,nbndf)
	        enddo  !! over ind_omega
	        deallocate(gx_omega)
		      call from_omega_to_tau_baa(gc_omega,nbndf,gc_tau,nbndf,n)
		    else if(way=='QP_based') then
		      call qp_green(ispin,ind_k)
		    endif
		    call output_gc_band_tau(ispin,ind_k,nbndf)
		    call output_gc_band_omega(ispin,k,nbndf)
		    call output_gx_band_tau(ispin,k,nbndf)
		    call output_sigc_band_tau(ispin,ind_k,nbndf)
	    endif
c ------ Contributions to E_tot --------------------------------------
	 	  call etot_gw_1(ispin,k,ind_k,gx_tau,n)
	    if(correlation) then  !! QP:GWG
	      g_full(:,:,:,:,ind_k,ispin)=gc_tau
	      do ind_tau=1,ndim3_tau
	        i_tau=me3_tau*ndim3_tau+ind_tau-1
	        ta=tau_mesh(i_tau)
              call g_x_tau(ispin,k,gx_tau,ta,n,chem_pot)
	        ta=betta_t-tau_mesh(i_tau)
              call g_x_tau(ispin,k,gx_tau(1,2),ta,n,chem_pot)
	        do i=1,n
	          cf(1)=gx_tau(i,1)-gx_tau(i,2)
	          cf(2)=gx_tau(i,1)+gx_tau(i,2)
	          do it=1,2
	            g_full(i,i,it,ind_tau,ind_k,ispin)=
     &	            g_full(i,i,it,ind_tau,ind_k,ispin)+cf(it)
	          enddo
	        enddo
	      enddo  !! over ind_tau
		  endif
	    deallocate(gx_tau)
          if(correlation) then
	      allocate(gx_omega(n))
	      do ind_omega=1,ndim3_omega
	        i_omega=me3_tau*ndim3_omega+ind_omega-1
		      call g_x_omega(ispin,k,gx_omega,w_omega(i_omega),n,
     &		                 chem_pot)
              call ferm_unpack_omega(tmp,gc_omega(1,1,1,ind_omega),
     &                               n,nbndf,nbndf)
              do i=1,n
                tmp(i,i)=tmp(i,i)+gx_omega(i)
              enddo
              call ferm_pack_omega(tmp,gc_omega(1,1,1,ind_omega),n,
     &                             nbndf,nbndf)
            enddo
            deallocate(gx_omega)
          endif
	  enddo  !! over ind_k
	enddo  !! over ispin
	deallocate(tmp)
	if(correlation) deallocate(gc_omega,gc_tau)
      call timel('******* G_KS_TAU finished **********')
c ------------- Forming G_FULL(0-) --> G_FULL_0 ------------------------
	call green_0(correlation)
c ------------- Forming Local Green's Function -------------------------
      call g_mt_from_g_full(correlation)
      end
