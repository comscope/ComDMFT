      subroutine ws_solid_resp(ind_k,w_last,w_red)
	use atom_mod
	use etot_mod
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	integer, intent(in) :: ind_k
	logical, intent(in) :: w_last,w_red
      integer :: k,n_pbt,i,ind_nu,i_nu,k0,info
	real*8 :: pi2a,q,q2,pi8,dd1,check_hermit
	complex*16 :: cc
	integer, allocatable :: ipiv(:)
	real*8, allocatable :: w_c_nu(:,:,:)
	complex*16, allocatable :: v_tmp(:,:),v_p(:,:),w_s(:,:),
     &                           work(:,:),work1(:,:)
      if(ind_k==1) call timel('***** Ws_solid_resp started ********')
	k=n3_mpi_k(me3_k+1)+ind_k
	k0=k_npnt_in_line(k)
	pi2a=(pi+pi)/par
	pi8=8.d0*pi
	q2=pi2a**2*dot_product(pnt(:,k),pnt(:,k))
	q=sqrt(q2)
	n_pbt=n_pbmt+nplwgw(k)
	allocate(ipiv(n_pbt))
	allocate(work(n_pbt,n_pbt))
	allocate(work1(n_pbt,n_pbt))
	allocate(w_s(n_pbt,n_pbt))
	allocate(v_tmp(n_pbt,n_pbt))
	call v_coul_full(nplwgw(k),v_tmp,ind_k)
c   ------- Transform R_FULL from TAU to NU representation for given Q -
  	allocate(w_c_nu(n_pbt,n_pbt,ndim3_nu))
	call from_tau_to_nu_ba(w_c_nu,n_pbt,p_w_solid_q_tau(1,1,1,ind_k),
     &	                     n_pbtot,n_pbt)
c --------------------------------------------------------------------
	allocate(v_p(n_pbt,n_pbt))
c --------------------------------------------------------------------
	dd1=0.d0
      if(ind_k==1) call timel('***** W: NU_loop is starting *******')
	do ind_nu=1,ndim3_nu
	  i_nu=me3_tau*ndim3_nu+ind_nu-1
	  call unpack_hermit(work,w_c_nu(1,1,ind_nu),n_pbt,n_pbt,n_pbt,
     &	                   (0.d0,0.d0),(1.d0,0.d0))
c --------------------------------------------------------------------
	  if(q2<1.d-10.and.v_coul_set==0) then
	    if(i_nu==0.and.metal) then
	      work1=work     !! R
	      call invers_h(n_pbt,work1,n_pbt)   !! R^(-1)
	      work1=work1+v_tmp                  !! R^(-1)+V
	      call invers_h(n_pbt,work1,n_pbt)   !! (R^(-1)+V)^(-1)
	      w_s=v_tmp
	      call invers_h(n_pbt,w_s,n_pbt)
	      w_s=w_s-work1
	      call invers_h(n_pbt,w_s,n_pbt)
	    else
            w_s=em1_gg0(1,i_nu)*v_tmp
          endif
c ------- We enforce the hermiticity here ------------------------------
	    w_s=w_s-v_tmp
          dd1=max(dd1,check_hermit(w_s,n_pbt,n_pbt))
	  else !!! if(k/=1) then
c  ------- V_P = V*R ------------------------------
	    call zgemm('n','n',n_pbt,n_pbt,n_pbt,(1.d0,0.d0),v_tmp,n_pbt,
     &	             work,n_pbt,(0.d0,0.d0),v_p,n_pbt)
c --------- V_P = 1 + V*R ------------------------------------------
	    do i=1,n_pbt
	      v_p(i,i)=(1.d0,0.d0)+v_p(i,i)
	    enddo
c --------- W_s = (1 + V*R) * V   -----------------------------
	    call zgemm('n','n',n_pbt,n_pbt,n_pbt,(1.d0,0.d0),v_p,n_pbt,
     &	             v_tmp,n_pbt,(0.d0,0.d0),w_s,n_pbt)
c ------- Contribution to eps^(-1) for q=(0,0,0) --------------------
          call pw_from_pb(ind_k,n_pbt,n_pbt,1,v_p,cc,2)
	    if(k0/=0) then
	      em1_gg0(k0,i_nu)=real(cc)
		    if(w_last) then
		      e1_output(i_nu,k0)=real(cc)
c -------- Store for output P ----------------------------------------
              work1=transpose(work)
              v_p=transpose(v_p)
              call zgesv(n_pbt,n_pbt,v_p,n_pbt,ipiv,work1,n_pbt,info)
              work1=transpose(work1)
              call pw_from_pb(ind_k,n_pbt,n_pbt,1,work1,cc,3)
	        p_output(i_nu,k0)=real(cc)
		    endif
		  endif
          w_s=w_s-v_tmp
	  endif
	  if(w_last) then
c ----- Contribution to Etot_C from PW --------------------------------	  
	    call etot_c_from_pw(n_pbt,k,i_nu,work,v_tmp)
c ----- Contribution to Phi_c in RPA  --------------------------------	  
	    call phi_rpa_c(n_pbt,n_pbt,ind_k,i_nu,work,v_tmp)
	  endif
c -------------------------------------------------------------------
        call pack_hermit(w_s,w_c_nu(1,1,ind_nu),n_pbt,n_pbt,n_pbt,0.d0,
     &                   1.d0)
        if(w_red) then
	    if(k==1.and..not.include_gamma_w) then
	      w_red_q(:,:,ind_nu,ind_k)=(0.d0,0.d0)
	    else
            call rpb_from_fpb_full(w_red_q(1,1,ind_nu,ind_k),w_s,k,
     &                             n_pbtot_red,n_pbt)
          endif
        endif
c ---------------------------------------------------------------------
	enddo  !! over ind_nu
	deallocate(v_p)
      if(q2<1.d-10) then
	  if(nproc_tau/=1) call dgop(dd1,1,'max',comm_pnk)
	  if(maswrk) then
	    write(iun,*)' Max deviation from Hermicity in W is ',dd1
	  endif
	endif
      if(ind_k==1) call timel('***** W: NU_loop finished **********')
c -------- Store for output W ----------------------------------------
      if(k0/=0) call add_to_output_w(ind_k,k0,n_pbt,w_c_nu,1)
c   ------- Transform NU ---> TAU representation ---------------------
      call from_nu_to_tau_ba(w_c_nu,n_pbt,
     &                       p_w_solid_q_tau(1,1,1,ind_k),n_pbtot,
     &                       n_pbt)
      if(ind_k==1) call timel('***** W: NU_to_TAU finished ********')
	deallocate(w_c_nu)
	deallocate(ipiv,work,w_s,work1,v_tmp)
      end
