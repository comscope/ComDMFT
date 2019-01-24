	subroutine transf_bos_from_tau(tran,nom,nom0,wom)
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
      integer, intent(in) :: nom,nom0
      real*8, intent(in) :: wom(nom)
      real*8, intent(out) :: tran(n_tau/2,nom0,2)
	integer :: k,int,l,n,n_int,i_nu,ind_int,nd,n1
	real*8 :: w_n,hh,x0,x1,ai,bi,aa,bb,yk,tk,xk,an,bn,ca,sa,cb,sb,cb0,
     &	        sb0,cb1,sb1,cb2,sb2,cb3,sb3,ssum3,csum2,ssum1,csum0
      integer, allocatable :: n_dim(:),n_mpi(:)
	real*8, allocatable :: s_bos(:,:),c_bos(:,:),cbar(:,:),tn(:),
     &	                     tl(:),cn(:),sn(:)
	allocate(n_dim(nproc_k))
	allocate(n_mpi(nproc_k))
c ---------------------------------------------------------------
	n_int=n_tau_int
	call size_shift_par(n_int,nproc_k,n_dim,n_mpi)
	hh=1.d0/n_int
c ---------------------------------------------------------------
      allocate(cbar(0:3,0:n_tau-2))
      allocate(tn(0:3))
      allocate(cn(0:3))
      allocate(sn(0:3))
      allocate(tl(0:n_tau-2))
	allocate(c_bos(0:n_tau/2-1,nom))
	c_bos=0.d0
	allocate(s_bos(0:n_tau/2-2,nom))
	s_bos=0.d0
      do ind_int=1,n_dim(me3_k+1)
	  int=n_mpi(me3_k+1)+ind_int
	  x1=int*hh-1.d0
	  x0=x1-hh
	  ai=0.5d0*betta_t*(1.d0+x0)**3
	  bi=0.5d0*betta_t*(1.d0+x1)**3
	  aa=0.5d0*(ai+bi)
	  bb=0.5d0*(bi-ai)
	  cbar=0.d0
	  do k=1,4
	    yk=cos(pi*(4.5d0-k)/4.d0)
	    call dcheb_recurr(yk,tn,3)
	    tk=aa+bb*yk
	    xk=(2.d0*tk/betta_t)**(1.d0/3.d0)-1.d0
	    call dcheb_recurr(xk,tl,n_tau-2)
	    do l=0,n_tau-2
	      do n=0,3
	        cbar(n,l)=cbar(n,l)+tl(l)*tn(n)
	      enddo
	    enddo
	  enddo
	  cbar(0,:)=cbar(0,:)/4.d0
	  cbar(1:3,:)=cbar(1:3,:)/2.d0
        do i_nu=1,nom
	    w_n=wom(i_nu)
	    an=aa*w_n
	    bn=bb*w_n
	    ca=cos(an)
	    sa=sin(an)
	    cb=cos(bn)
	    sb=sin(bn)
c ------------------------------------------------------------------
	    sb0=0.d0
	    cb1=0.d0
	    sb2=0.d0
	    cb3=0.d0
	    if(abs(bn)<1.d-1) then
	      cb0=2.d0*csum0(bn)
	      sb1=2.d0*bn*ssum1(bn)
	      cb2=-4.d0*csum2(bn)-2.d0/3.d0
	      sb3=4.d0*ssum3(bn)
	    else	    
	      cb0=2.d0*sb/bn
	      sb1=2.d0*(sb-bn*cb)/bn**2
	      cb2=2.d0*((bn*bn-4.d0)*sb+4.d0*bn*cb)/bn**3
	      sb3=2.d0*((9.d0*bn*bn-24.d0)*sb-(bn**3-24.d0*bn)*cb)/bn**4
	    endif
c ------------------------------------------------------------------
	    cn(0)=ca*cb0-sa*sb0
	    cn(1)=ca*cb1-sa*sb1
	    cn(2)=ca*cb2-sa*sb2
	    cn(3)=ca*cb3-sa*sb3
	    sn(0)=sa*cb0+ca*sb0
	    sn(1)=sa*cb1+ca*sb1
	    sn(2)=sa*cb2+ca*sb2
	    sn(3)=sa*cb3+ca*sb3
	    do n=0,3
	      do l=0,n_tau/2-1
	        c_bos(l,i_nu)=c_bos(l,i_nu)+bb*cbar(n,l+l)*cn(n)
	      enddo
	      do l=0,n_tau/2-2
	        s_bos(l,i_nu)=s_bos(l,i_nu)+bb*cbar(n,l+l+1)*sn(n)
	      enddo
	    enddo
	  enddo   !! over i_nu
	enddo   !! over ind_int
      if(nproc_k/=1) then
        call DGOP(c_bos,n_tau/2*nom,'  +',comm_pnt)
        call DGOP(s_bos,(n_tau/2-1)*nom,'  +',comm_pnt)
      endif
c -------------------------------------------------------------
	call size_shift_par(nom,nproc_k,n_dim,n_mpi)
	tran=0.d0
	nd=n_dim(me3_k+1)
	n1=n_mpi(me3_k+1)+1
	if(nd/=0) then
	  call dgemm('t','n',n_tau/2,nd,n_tau/2,1.d0,
     &	           cheb_tau_even(0,1),n_tau/2,c_bos(0,n1),n_tau/2,0.d0,
     &             tran(1,n1,1),n_tau/2)
	  call dgemm('t','n',n_tau/2-1,nd,n_tau/2-1,1.d0,
     &	           cheb_tau_odd(0,1),n_tau/2-1,s_bos(0,n1),n_tau/2-1,
     &             0.d0,tran(1,n1,2),n_tau/2)
      endif
      if(nproc_k/=1) call DGOP(tran,n_tau*nom0,'  +',comm_pnt)
	deallocate(n_dim,n_mpi)
	deallocate(s_bos,c_bos,cbar,tn,tl,cn,sn)
      end