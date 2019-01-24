      subroutine transf_tau_from_nu_cheb
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
	integer :: ind_tau,i_nu,i_tau,n,ic_nu
	real*8 :: c,c1,cs,sn,w0,wm,a,b,dd,w,x,w1,c11,c12,c21,c22,s11,s12,
     &	        s21,s22
	real*8, allocatable :: cheb(:),tmp1(:,:,:)
c ----------------------------------------------------------------
      w0=sqrt(w_nu(n_nu_exa))
      wm=sqrt(nu_max)
      a=2.d0*w0*wm
      b=w0+wm
      c=wm-w0
	allocate(cheb(0:n_nu_asy-2))
	allocate(tmp1(0:n_nu_asy-2,ndim3_tau,2))
	tmp1=0.d0
	dd=2.d0*pi/betta_t
	w=w_nu(n_nu_exa)
	do i_nu=n_nu_exa+1,n_nu_maz
	  w=w+dd
        x=(b-a/sqrt(w))/c
	  call dcheb_recurr(x,cheb,n_nu_asy-2)
	  do ind_tau=1,ndim3_tau
	    i_tau=me3_tau*ndim3_tau+ind_tau-1
	    c1=w*tau_mesh(i_tau)
	    cs=cos(c1)
	    sn=sin(c1)
	    do n=0,n_nu_asy-2
	      tmp1(n,ind_tau,1)=tmp1(n,ind_tau,1)+cs*cheb(n)
	      tmp1(n,ind_tau,2)=tmp1(n,ind_tau,2)+sn*cheb(n)
	    enddo
	  enddo
	enddo
	deallocate(cheb)
      tau_from_nu=0.d0
      do i_nu=0,n_nu
	  w=w_nu(i_nu)
        if(i_nu==0) then
          tau_from_nu(i_nu,:,1)=0.5d0
          tau_from_nu(i_nu,:,2)=0.d0
        else if(i_nu<=n_nu_exa) then
	    do ind_tau=1,ndim3_tau
	      i_tau=me3_tau*ndim3_tau+ind_tau-1
	      c1=w*tau_mesh(i_tau)
	      tau_from_nu(i_nu,ind_tau,1)=cos(c1)
	      tau_from_nu(i_nu,ind_tau,2)=sin(c1)
	    enddo
        else if(i_nu<=n_nu-1) then
          ic_nu=i_nu-n_nu_exa
	    do ind_tau=1,ndim3_tau
	      tau_from_nu(i_nu,ind_tau,1)=dot_product(tmp1(:,ind_tau,1),
     &	                                            cheb_nu(:,ic_nu))
	      tau_from_nu(i_nu,ind_tau,2)=dot_product(tmp1(:,ind_tau,2),
     &	                                            cheb_nu(:,ic_nu))
	    enddo
	  endif
	enddo
      deallocate(tmp1)
c ---- Here are Nu-Sums for Internal Energy evaluation ---------
      allocate(et_sum_nu(0:n_nu))
	if(me3_tau==0) et_sum_nu=tau_from_nu(:,1,1)
	if(nproc_tau/=1) call brdcst(msgdbl,et_sum_nu,8*(n_nu+1),0,
     &	                           comm_pnk)
c ------- Asymptotic contribution for NU ---> TAU --------------------
      w=w_nu(n_nu)
      w1=w_nu(n_nu-1)
      c11=w1**4/(w1-w)/(w1+w)
      c21=w**4/(w-w1)/(w1+w)
      c12=-w1**4*w**2/(w1-w)/(w1+w)
      c22=w1**2*w**4/(w1-w)/(w1+w)
      s11=w1**5/(w1-w)/(w1+w)
      s21=w**5/(w-w1)/(w1+w)
      s12=-w1**5*w**2/(w1-w)/(w1+w)
      s22=w1**2*w**5/(w1-w)/(w1+w)
	do ind_tau=1,ndim3_tau
	  i_tau=me3_tau*ndim3_tau+ind_tau-1
        tau_from_nu(n_nu-1,ind_tau,1)=tau_from_nu(n_nu-1,ind_tau,1)
     &                               +c11*sum_nu_asy(i_tau,1,1)
     &                               +c12*sum_nu_asy(i_tau,2,1)
        tau_from_nu(n_nu,ind_tau,1)=tau_from_nu(n_nu,ind_tau,1)
     &                             +c21*sum_nu_asy(i_tau,1,1)
     &                             +c22*sum_nu_asy(i_tau,2,1)
        tau_from_nu(n_nu-1,ind_tau,2)=tau_from_nu(n_nu-1,ind_tau,2)
     &                               +s11*sum_nu_asy(i_tau,2,2)
     &                               +s12*sum_nu_asy(i_tau,3,2)
        tau_from_nu(n_nu,ind_tau,2)=tau_from_nu(n_nu,ind_tau,2)
     &                             +s21*sum_nu_asy(i_tau,2,2)
     &                             +s22*sum_nu_asy(i_tau,3,2)
      enddo
c ------- Asymptotic contribution for Etotal = P*W --------------------
      c11=w1**6/(w1-w)/(w1+w)
      c21=w**6/(w-w1)/(w1+w)
      c12=-w1**6*w**2/(w1-w)/(w1+w)
      c22=w1**2*w**6/(w1-w)/(w1+w)
      et_sum_nu(n_nu-1)=et_sum_nu(n_nu-1)+c11*sum_nu_asy(0,2,1)
     &                                   +c12*sum_nu_asy(0,3,1)
      et_sum_nu(n_nu)=et_sum_nu(n_nu)+c21*sum_nu_asy(0,2,1)
     &                               +c22*sum_nu_asy(0,3,1)
c ---------- Normalization -----------------------------------------
      tau_from_nu=tau_from_nu*4.d0/betta_t
      et_sum_nu=et_sum_nu/betta_t
      end