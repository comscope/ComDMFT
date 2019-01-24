	subroutine transf_bos_from_tau1(tran,n,n0,i_tau,cheb_tau,key)
c     key = 1 integral from 0 to tau
c     key = 2 integral from tau to beta	
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
      integer, intent(in) :: n,n0,i_tau,key
      real*8, intent(in) :: cheb_tau(0:n_tau-2,n_tau-1,2)
      real*8, intent(out) :: tran(0:n_tau,0:n_nu_d,2,2)
	integer :: k,int,l,n_int,i_nu,m_tau,i
	real*8 :: w_n,hh,x0,x1,ai,bi,aa,bb,yk,tk,xk,an,bn,ca,sa,cb,sb,cb0,
     &	        sb0,cb1,sb1,cb2,sb2,cb3,sb3,a,ssum3,csum2,ssum1,csum0
	real*8, allocatable :: s_bos(:,:),c_bos(:,:),cbar(:,:),tn(:),
     &	                     tl(:),cn(:),sn(:)
c ---------------------------------------------------------------
	n_int=n_tau_int
	hh=1.d0/n_int
        m_tau=0
	if(key==1) m_tau=n
	if(key==2) m_tau=n0-n
	if(m_tau==0) return
c ---------------------------------------------------------------
      allocate(cbar(0:3,0:m_tau-2))
      allocate(tn(0:3))
      allocate(cn(0:3))
      allocate(sn(0:3))
      allocate(tl(0:m_tau-2))
	allocate(c_bos(0:m_tau-2,0:n_nu_d))
	c_bos=0.d0
	allocate(s_bos(0:m_tau-2,0:n_nu_d))
	s_bos=0.d0
      do int=1,n_int
	  x1=int*hh-1.d0
	  x0=x1-hh
	  if(key==1) then
	    ai=0.5d0*tau_mesh(i_tau)*(1.d0+x0)**3
	    bi=0.5d0*tau_mesh(i_tau)*(1.d0+x1)**3
	  else
          a=0.5d0*(betta_t-tau_mesh(i_tau))
	    ai=tau_mesh(i_tau)+a*(1.d0+x0)**3
	    bi=tau_mesh(i_tau)+a*(1.d0+x1)**3
	  endif
	  aa=0.5d0*(ai+bi)
	  bb=0.5d0*(bi-ai)
	  cbar=0.d0
	  do k=1,4
	    yk=cos(pi*(4.5d0-k)/4.d0)
	    call dcheb_recurr(yk,tn,3)
	    tk=aa+bb*yk
	    if(key==1) then
	      xk=(2.d0*tk/tau_mesh(i_tau))**(1.d0/3.d0)-1.d0
	    else if(key==2) then
	      xk=(2.d0*(tk-tau_mesh(i_tau))
     &	      /(betta_t-tau_mesh(i_tau)))**(1.d0/3.d0)-1.d0
	    endif
	    call dcheb_recurr(xk,tl,m_tau-2)
	    do l=0,m_tau-2
	      do i=0,3
	        cbar(i,l)=cbar(i,l)+tl(l)*tn(i)
	      enddo
	    enddo
	  enddo
	  cbar(0,:)=cbar(0,:)/4.d0
	  cbar(1:3,:)=cbar(1:3,:)/2.d0
        do i_nu=0,n_nu_d
	    w_n=w_nu_d(i_nu)
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
	    do i=0,3
	      do l=0,m_tau-2
	        c_bos(l,i_nu)=c_bos(l,i_nu)+bb*cbar(i,l)*cn(i)
	        s_bos(l,i_nu)=s_bos(l,i_nu)+bb*cbar(i,l)*sn(i)
	      enddo
	    enddo
	  enddo   !! over i_nu
	enddo   !! over int
c -------------------------------------------------------------
      if(key==1) then
c ----------------- A1 COS ------------------------------------
	  call dgemm('t','n',m_tau-1,n_nu_d+1,m_tau-1,1.d0,
     & 	           cheb_tau(0,1,1),n_tau-1,c_bos(0,0),m_tau-1,0.d0,
     &             tran(1,0,1,1),n_tau+1)
c ----------------- A1 SIN ------------------------------------
	  call dgemm('t','n',m_tau-1,n_nu_d+1,m_tau-1,1.d0,
     & 	           cheb_tau(0,1,1),n_tau-1,s_bos(0,0),m_tau-1,0.d0,
     &             tran(1,0,2,1),n_tau+1)
c ----------------- B1 COS ------------------------------------
	  call dgemm('t','n',m_tau-1,n_nu_d+1,m_tau-1,1.d0,
     & 	           cheb_tau(0,1,2),n_tau-1,c_bos(0,0),m_tau-1,0.d0,
     &             tran(1,0,1,2),n_tau+1)
c ----------------- B1 SIN ------------------------------------
	  call dgemm('t','n',m_tau-1,n_nu_d+1,m_tau-1,1.d0,
     & 	           cheb_tau(0,1,2),n_tau-1,s_bos(0,0),m_tau-1,0.d0,
     &             tran(1,0,2,2),n_tau+1)
      else if(key==2) then
c ----------------- A2 COS ------------------------------------
	  call dgemm('t','n',m_tau-1,n_nu_d+1,m_tau-1,1.d0,
     & 	           cheb_tau(0,1,1),n_tau-1,c_bos(0,0),m_tau-1,0.d0,
     &             tran(n+1,0,1,1),n_tau+1)
c ----------------- A2 SIN ------------------------------------
	  call dgemm('t','n',m_tau-1,n_nu_d+1,m_tau-1,1.d0,
     & 	           cheb_tau(0,1,1),n_tau-1,s_bos(0,0),m_tau-1,0.d0,
     &             tran(n+1,0,2,1),n_tau+1)
c ----------------- B2 COS ------------------------------------
	  call dgemm('t','n',m_tau-1,n_nu_d+1,m_tau-1,1.d0,
     & 	           cheb_tau(0,1,2),n_tau-1,c_bos(0,0),m_tau-1,0.d0,
     &             tran(n+1,0,1,2),n_tau+1)
c ----------------- B2 SIN ------------------------------------
	  call dgemm('t','n',m_tau-1,n_nu_d+1,m_tau-1,1.d0,
     & 	           cheb_tau(0,1,2),n_tau-1,s_bos(0,0),m_tau-1,0.d0,
     &             tran(n+1,0,2,2),n_tau+1)
      endif
	deallocate(s_bos,c_bos,cbar,tn,tl,cn,sn)
      end
