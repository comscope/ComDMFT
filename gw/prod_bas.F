      subroutine prod_bas
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: limax,ndim0,mx2,nli,k
	integer, allocatable :: nl(:,:),ind_l(:,:,:)
	limax=nrel*maxb+1
c	key=0 - Run to find the dimensions
	allocate(nl(limax,limax))
	allocate(ind_l(limax,limax,maxpb+1))
	call prod_b1(ndim0,nl,ind_l,limax)
      if(maswrk) write(iun,*)' Max dimension in full PB search is ',
     &                       ndim0
	  call prod_b2(ndim0,0,nl,ind_l,limax)
	  if(allocated(ind_pb)) deallocate(ind_pb)
	  allocate(ind_pb(maxlfpb,nsort))
	  if(allocated(lval_pb)) deallocate(lval_pb)
	  allocate(lval_pb(maxlfpb,nsort))
	  if(allocated(fun_pb)) deallocate(fun_pb)
	  allocate(fun_pb(maxpbr))
	  if(allocated(dfun_pb)) deallocate(dfun_pb)
	  allocate(dfun_pb(max_e_pb,nsort))
	  if(allocated(lm_pbmt)) deallocate(lm_pbmt)
	  allocate(lm_pbmt(n_pbmtm,nsort))
	  if(allocated(lf_pb)) deallocate(lf_pb)
	  allocate(lf_pb(n_pbmtm,nsort))
	  if(allocated(ind_prod0)) deallocate(ind_prod0)
	  allocate(ind_prod0(max_e_pb,0:maxpb,nsort))
	  if(allocated(ind_prod)) deallocate(ind_prod)
	  allocate(ind_prod(max_e_pb,limlpb,nsort))
c	  key=1 - Run to fill the arrays
	  call prod_b2(ndim0,1,nl,ind_l,limax)
	  if(allocated(fi0)) deallocate(fi0)
        allocate(fi0(maxlfun,maxlfun,maxlfpb,nsort,nspin,nspin))
	  if(allocated(tmm_pb)) deallocate(tmm_pb)
	  allocate(tmm_pb(max_e_pb,max_e_pb,0:maxpb,nsort))
	  call get_fi0_gw
	  call get_tmm
	  n_pbtot=n_pbmt+nplw_gw
	  if(.not.allocated(ndim_pbmt)) allocate(ndim_pbmt(nproc_k))
	  if(.not.allocated(n_mpi_pbmt)) allocate(n_mpi_pbmt(nproc_k))
	  call size_shift_par(n_pbmt,nproc_k,ndim_pbmt,n_mpi_pbmt)
	  if(.not.allocated(ndim_pbtt)) allocate(ndim_pbtt(nproc_tau))
	  if(.not.allocated(n_mpi_pbtt)) allocate(n_mpi_pbtt(nproc_tau))
	  call size_shift_par(n_pbmt,nproc_tau,ndim_pbtt,n_mpi_pbtt)
	  if(.not.allocated(ndim_pbt)) allocate(ndim_pbt(nproc_tau))
	  if(.not.allocated(n_mpi_pbt)) allocate(n_mpi_pbt(nproc_tau))
	  call size_shift_par(n_pbtot,nproc_tau,ndim_pbt,n_mpi_pbt)
	  if(maswrk) then
	    write(iun,*)' The size of PB is (k-point list) :'
	    write(iun,'(14i5)')(n_pbmt+nplwgw(k),k=1,npnt)
	    write(iun,'(a22,i6,a1,i3,a1,a5,i6,a7,i6)')
     &	               'Product basis: MT-part',n_pbmt,'(',nlb2,')',
     &                   ' Int:',nplw_gw,' Total:',n_pbtot
	  endif
c ------------- Coulomb Interaction matrix -------------------------
	  if(allocated(v_coul_pb)) deallocate(v_coul_pb)
	  allocate(v_coul_pb(n_pbtot,ndim_pbt(me3_tau+1),
     &             maxval(ndim3_k)))
        call v_coulomb 
	  call v_im_r_from_k(0)    
	deallocate(nl,ind_l)
      call timel('**** prod_bas finished *************')
      end
