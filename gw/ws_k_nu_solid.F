      subroutine ws_k_nu_solid(w_lda)
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	logical, intent(in) :: w_lda
	logical :: w_last,w_red
      integer :: k,ind_k,npbt1,mx2,it0
	real*8 :: ak,dev
	mx2=maxel_loc**2
	w_last=.false.
	it0=iter_dft+iter_hf+iter_gw
	if(ubi==' gw'.and.iter==it0) w_last=.true.
	if(ubi==' qp'.and.iter==it0+iter_qp) w_last=.true.
      w_red=.false.
	if(w_last) then
	  if(allocated(e_output)) deallocate(e_output)
	  allocate(e_output(0:n_nu,n_line))
	  e_output=0.d0
	  allocate(e1_output(0:n_nu,n_line))
	  e1_output=0.d0
	endif
	allocate(w_output(0:n_nu,n_line))
	w_output=0.d0
	allocate(p_output(0:n_nu,n_line))
	p_output=0.d0
	allocate(em1_gg0(n_line,0:n_nu))
	em1_gg0=0.d0
	npbt1=n_pbmt+nplwgw(1)
c -----------------------------------------------------------------------
	dev=0.d0
	e_c_solid_pw=0.d0
	phi_functional_c=0.d0
	do ind_k=1,ndim3_k(me3_k+1)
	  k=n3_mpi_k(me3_k+1)+ind_k
	  ak=dot_product(pnt(:,k),pnt(:,k))
	  if(ak>1.d-12) then
	      call ws_k_nu_solid_0(ind_k,w_last,w_red)
	  endif
	enddo  !! over k
	if(v_coul_set==0) then
c --------- Collect EM1_gg0 ------------------------------------------
	  if(nproc_k/=1) then
	    call dgop(em1_gg0,n_line*(n_nu+1),'  +',comm_pnt)
	  endif
	  if(nproc_tau/=1) call dgop(em1_gg0,n_line*(n_nu+1),'  +',
     &	                           comm_pnk)
	  if(.not.metal.and.n_line<=5) call e1_q0_fit(em1_gg0)
	  if(.not.metal.and.n_line>5) call e2_q0_fit(em1_gg0)
	  if(metal) call e2_q0_fit(em1_gg0)
	endif
c --------------------------------------------------------------------
	if(me3_k==0) then
	    call ws_k_nu_solid_0(1,w_last,w_red)
	endif
	if(nproc_k/=1) then
	  call dgop(e_c_solid_pw,1,'  +',comm_pnt)
	  call dgop(phi_functional_c,1,'  +',comm_pnt)
	  call dgop(w_output,n_line*(n_nu+1),'  +',comm_pnt)
	  call dgop(p_output,n_line*(n_nu+1),'  +',comm_pnt)
	  if(w_last) then
	    call dgop(e_output,n_line*(n_nu+1),'  +',comm_pnt)
	    call dgop(e1_output,n_line*(n_nu+1),'  +',comm_pnt)
	  endif
	endif
	if(nproc_tau/=1) then
	  call dgop(e_c_solid_pw,1,'  +',comm_pnk)
	  call dgop(phi_functional_c,1,'  +',comm_pnk)
	  call dgop(w_output,n_line*(n_nu+1),'  +',comm_pnk)
	  call dgop(p_output,n_line*(n_nu+1),'  +',comm_pnk)
	  if(w_last) then
	    call dgop(e_output,n_line*(n_nu+1),'  +',comm_pnk)
	    call dgop(e1_output,n_line*(n_nu+1),'  +',comm_pnk)
	  endif
	endif
	if(v_coul_set==0.and.w_last) e1_output(:,1)=em1_gg0(1,:)
c -------------------------------------------------------------- 
	call output_w(w_last)
	deallocate(em1_gg0)
	deallocate(w_output,p_output)
	if(w_last) then
	  deallocate(e_output,e1_output)
      endif
      end
