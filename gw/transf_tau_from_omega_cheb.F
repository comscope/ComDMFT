      subroutine transf_tau_from_omega_cheb
c     Creates transformation matrices for OMEGA ---> TAU
c       They correspond to the transforms :
c     A(t) = 4/b Sum_w>0 cos(wt) A(w)
c     B(t) = 4/b Sum_w>0 sin(wt) B(w)
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
	integer :: ind_tau,i_omega,i_tau,n,ic_omega
	real*8 :: c,c1,cs,sn,w0,wm,a,b,dd,w,x,w1,c11,c12,c21,c22,s11,s12,
     &	        s21,s22
	real*8, allocatable :: cheb(:),tmp1(:,:,:)
c ----------------------------------------------------------------
      w0=sqrt(w_omega(n_omega_exa))
      wm=sqrt(omega_max)
      a=2.d0*w0*wm
      b=w0+wm
      c=wm-w0
	allocate(cheb(0:n_omega_asy-2))
	allocate(tmp1(0:n_omega_asy-2,ndim3_tau,2))
	tmp1=0.d0
	dd=2.d0*pi/betta_t
	w=w_omega(n_omega_exa)
	do i_omega=n_omega_exa+1,n_omega_maz
	  w=w+dd
        x=(b-a/sqrt(w))/c
	  call dcheb_recurr(x,cheb,n_omega_asy-2)
	  do ind_tau=1,ndim3_tau
	    i_tau=me3_tau*ndim3_tau+ind_tau-1
	    c1=w*tau_mesh(i_tau)
	    cs=cos(c1)
	    sn=sin(c1)
	    do n=0,n_omega_asy-2
	      tmp1(n,ind_tau,1)=tmp1(n,ind_tau,1)+cs*cheb(n)
	      tmp1(n,ind_tau,2)=tmp1(n,ind_tau,2)+sn*cheb(n)
	    enddo
	  enddo
	enddo
	deallocate(cheb)
      tau_from_omega=0.d0
      do i_omega=0,n_omega
	  w=w_omega(i_omega)
        if(i_omega<=n_omega_exa) then
	    do ind_tau=1,ndim3_tau
	      i_tau=me3_tau*ndim3_tau+ind_tau-1
	      c1=w*tau_mesh(i_tau)
	      tau_from_omega(i_omega,ind_tau,1)=cos(c1)
	      tau_from_omega(i_omega,ind_tau,2)=sin(c1)
	    enddo
        else if(i_omega<=n_omega-1) then
          ic_omega=i_omega-n_omega_exa
	    do ind_tau=1,ndim3_tau
	      tau_from_omega(i_omega,ind_tau,1)=
     &	      dot_product(tmp1(:,ind_tau,1),cheb_omega(:,ic_omega))
	      tau_from_omega(i_omega,ind_tau,2)=
     &	      dot_product(tmp1(:,ind_tau,2),cheb_omega(:,ic_omega))
	    enddo
	  endif
	enddo
      deallocate(tmp1)
c ---- Here are Omega-Sums for Internal Energy evaluation ---------
      allocate(et_sum_sc_gx(0:n_omega))
      allocate(et_sum_sc_gc(0:n_omega))
	if(me3_tau==0) et_sum_sc_gx=tau_from_omega(:,1,1)
	if(nproc_tau/=1) call brdcst(msgdbl,et_sum_sc_gx,8*(n_omega+1),0,
     &	                           comm_pnk)
      et_sum_sc_gc=et_sum_sc_gx  !! they differ only in asymptotic
c ------- Asymptotic contribution for OMEGA ---> TAU -----------------
      w=w_omega(n_omega)
      w1=w_omega(n_omega-1)
      c11=w1**6/(w1-w)/(w1+w)
      c21=w**6/(w-w1)/(w1+w)
      c12=-w1**6*w**2/(w1-w)/(w1+w)
      c22=w1**2*w**6/(w1-w)/(w1+w)
      s11=w1**5/(w1-w)/(w1+w)
      s21=w**5/(w-w1)/(w1+w)
      s12=-w1**5*w**2/(w1-w)/(w1+w)
      s22=w1**2*w**5/(w1-w)/(w1+w)
	do ind_tau=1,ndim3_tau
	  i_tau=me3_tau*ndim3_tau+ind_tau-1
        tau_from_omega(n_omega-1,ind_tau,1)=
     &    tau_from_omega(n_omega-1,ind_tau,1)
     &   +c11*sum_omega_asy(i_tau,2,1)
     &   +c12*sum_omega_asy(i_tau,3,1)
        tau_from_omega(n_omega,ind_tau,1)=
     &    tau_from_omega(n_omega,ind_tau,1)
     &   +c21*sum_omega_asy(i_tau,2,1)
     &   +c22*sum_omega_asy(i_tau,3,1)
        tau_from_omega(n_omega-1,ind_tau,2)=
     &    tau_from_omega(n_omega-1,ind_tau,2)
     &   +s11*sum_omega_asy(i_tau,2,2)
     &   +s12*sum_omega_asy(i_tau,3,2)
        tau_from_omega(n_omega,ind_tau,2)=
     &    tau_from_omega(n_omega,ind_tau,2)
     &   +s21*sum_omega_asy(i_tau,2,2)
     &   +s22*sum_omega_asy(i_tau,3,2)
      enddo  !! over ind_tau
c ------- Asymptotic contribution for Etotal = Sigma_C*G_C -------------
      c11=w1**6/(w1-w)/(w1+w)
      c21=w**6/(w-w1)/(w1+w)
      c12=-w1**6*w**2/(w1-w)/(w1+w)
      c22=w1**2*w**6/(w1-w)/(w1+w)
      et_sum_sc_gc(n_omega-1)=et_sum_sc_gc(n_omega-1)
     &                       +c11*sum_omega_asy(0,2,1)
     &                       +c12*sum_omega_asy(0,3,1)
      et_sum_sc_gc(n_omega)=et_sum_sc_gc(n_omega)
     &                     +c21*sum_omega_asy(0,2,1)
     &                     +c22*sum_omega_asy(0,3,1)
c ------- Asymptotic contribution for Etotal = Sigma_C*G_X -------------
      c11=w1**4/(w1-w)/(w1+w)
      c21=w**4/(w-w1)/(w1+w)
      c12=-w1**4*w**2/(w1-w)/(w1+w)
      c22=w1**2*w**4/(w1-w)/(w1+w)
      et_sum_sc_gx(n_omega-1)=et_sum_sc_gx(n_omega-1)
     &                       +c11*sum_omega_asy(0,1,1)
     &                       +c12*sum_omega_asy(0,2,1)
      et_sum_sc_gx(n_omega)=et_sum_sc_gx(n_omega)
     &                     +c21*sum_omega_asy(0,1,1)
     &                     +c22*sum_omega_asy(0,2,1)
c ---------- Normalization -----------------------------------------
      tau_from_omega=tau_from_omega*4.d0/betta_t
      et_sum_sc_gc=et_sum_sc_gc/betta_t
      et_sum_sc_gx=et_sum_sc_gx/betta_t
      end