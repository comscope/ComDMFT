      subroutine ferm_sums
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
	integer :: ind_tau,i_omega,i,i_tau
	real*8 :: fr,c1,cs,sn,tmp,aa,t,w,w2,w3,w4,w5,w6,
     &	        pi_bet,h1,nu1,nu2,nu3,dt,a1,a2,a3,b1,b2,b3,c2,c3,x4,x6,
     &          x8,v,v2,v4,v6,v8,ss2,ss4,ss6,ss8
	pi_bet=pi/betta_t
	h1=pi_bet+pi_bet
c ---- Here are Omega-Sums for Internal Energy evaluation ---------
      allocate(et_sum_sc_gc(0:n_omega))
      allocate(et_sum_sc_gx(0:n_omega))
      et_sum_sc_gc=0.d0
      et_sum_sc_gx=0.d0
c ---------------------- Exact Part ------------------------------
	do i_omega=0,n_omega_exa
	  et_sum_sc_gc(i_omega)=1.d0
	  et_sum_sc_gx(i_omega)=1.d0
	enddo
c ---- Auxiliary sums --------------------------------------------
      ss2=0.d0
      ss4=0.d0
      ss6=0.d0
      ss8=0.d0
      do i_omega=1001,mats_num_omega(n_omega-2)
        w=(2*i_omega+1)*pi/betta_t
        w2=w*w
        ss2=ss2+1.d0/w2	
        ss4=ss4+1.d0/w2**2	
        ss6=ss6+1.d0/w2**3
        ss8=ss8+1.d0/w2**4
      enddo	
c ----------------- Asymptotic Part ------------------------------	  
	do i=n_omega_exa+1,n_omega-1,2
	  nu1=w_omega(i-1)
	  nu2=w_omega(i)
	  nu3=w_omega(i+1)
c ---------------------- Transformation for G_c -----------------------	    
	  dt=((1.d0/nu3**2-1.d0/nu2**2)/(nu2*nu3)**2
     &     +(1.d0/nu2**2-1.d0/nu1**2)/(nu1*nu2)**2
     &	   +(1.d0/nu1**2-1.d0/nu3**2)/(nu1*nu3)**2)/(nu1*nu2*nu3)**4
	  a1=(1.d0/nu3**2-1.d0/nu2**2)/(nu2*nu3)**6/dt
	  a2=(1.d0/nu1**2-1.d0/nu3**2)/(nu1*nu3)**6/dt
	  a3=(1.d0/nu2**2-1.d0/nu1**2)/(nu1*nu2)**6/dt
	  b1=(1.d0/nu2**4-1.d0/nu3**4)/(nu2*nu3)**4/dt
	  b2=(1.d0/nu3**4-1.d0/nu1**4)/(nu1*nu3)**4/dt
	  b3=(1.d0/nu1**4-1.d0/nu2**4)/(nu1*nu2)**4/dt
	  c1=(1.d0/nu3**2-1.d0/nu2**2)/(nu2*nu3)**4/dt
	  c2=(1.d0/nu1**2-1.d0/nu3**2)/(nu1*nu3)**4/dt
	  c3=(1.d0/nu2**2-1.d0/nu1**2)/(nu1*nu2)**4/dt
	  x4=0.d0
	  x6=0.d0
	  x8=0.d0
	  do i_omega=mats_num_omega(i-1)+1,mats_num_omega(i+1)
	    v=h1*i_omega
	    v2=v*v
	    v4=v2*v2
	    v6=v4*v2
	    v8=v4*v4
	    x4=x4+1.d0/v4
	    x6=x6+1.d0/v6
	    x8=x8+1.d0/v8
	  enddo
	  et_sum_sc_gc(i-1)=et_sum_sc_gc(i-1)+x4*a1+x6*b1+x8*c1
	  et_sum_sc_gc(i)=et_sum_sc_gc(i)+x4*a2+x6*b2+x8*c2
	  et_sum_sc_gc(i+1)=et_sum_sc_gc(i+1)+x4*a3+x6*b3+x8*c3
c ---------------------- Transformation for G_x -----------------------	    
	  dt=((1.d0/nu3**2-1.d0/nu2**2)/(nu2*nu3)**2
     &     +(1.d0/nu2**2-1.d0/nu1**2)/(nu1*nu2)**2
     &	   +(1.d0/nu1**2-1.d0/nu3**2)/(nu1*nu3)**2)/(nu1*nu2*nu3)**2
	  a1=(1.d0/nu3**2-1.d0/nu2**2)/(nu2*nu3)**4/dt
	  a2=(1.d0/nu1**2-1.d0/nu3**2)/(nu1*nu3)**4/dt
	  a3=(1.d0/nu2**2-1.d0/nu1**2)/(nu1*nu2)**4/dt
	  b1=(1.d0/nu2**4-1.d0/nu3**4)/(nu2*nu3)**2/dt
	  b2=(1.d0/nu3**4-1.d0/nu1**4)/(nu1*nu3)**2/dt
	  b3=(1.d0/nu1**4-1.d0/nu2**4)/(nu1*nu2)**2/dt
	  c1=(1.d0/nu3**2-1.d0/nu2**2)/(nu2*nu3)**2/dt
	  c2=(1.d0/nu1**2-1.d0/nu3**2)/(nu1*nu3)**2/dt
	  c3=(1.d0/nu2**2-1.d0/nu1**2)/(nu1*nu2)**2/dt
	  x4=0.d0
	  x6=0.d0
	  x8=0.d0
	  if(i<n_omega-1) then
	    do i_omega=mats_num_omega(i-1)+1,mats_num_omega(i+1)
	      v=h1*i_omega
	      v2=v*v
	      v4=v2*v2
	      v6=v4*v2
	      x4=x4+1.d0/v2
	      x6=x6+1.d0/v4
	      x8=x8+1.d0/v6
	    enddo
	  else
c	    x4=betta_t**2/8.d0-ss2
c	    x6=betta_t**4/96.d0-ss4
c	    x8=betta_t**6/960.d0-ss6
	    x4=betta_t**2*2.49750229d-4/pi**2-ss2
	    x6=betta_t**4*2.077094776058d-11/pi**4-ss4
	    x8=betta_t**6*3.109417886834d-18/pi**6-ss6
	  endif
	  et_sum_sc_gx(i-1)=et_sum_sc_gx(i-1)+x4*a1+x6*b1+x8*c1
	  et_sum_sc_gx(i)=et_sum_sc_gx(i)+x4*a2+x6*b2+x8*c2
	  et_sum_sc_gx(i+1)=et_sum_sc_gx(i+1)+x4*a3+x6*b3+x8*c3
	enddo
	et_sum_sc_gc=et_sum_sc_gc/betta_t	
	et_sum_sc_gx=et_sum_sc_gx/betta_t	
c
c ------- Here sums for asymptotic contribution to TAU-functions ---
c
      sums_omega_asymp=0.d0
	do ind_tau=1,ndim_tau(me+1)
	  i_tau=n_mpi_tau(me+1)+ind_tau-1
        t=tau_mesh(i_tau)
        fr=pi/betta_t
c ---- We use some known sums here for cos/w; sin/w; cos/w2; ...----       
        do i=0,n_omega_maz
          w=(i+i+1.d0)*fr
          w2=w*w
          w3=w2*w
          w4=w2*w2
          w5=w3*w2
          tmp=w*t
          cs=cos(tmp)
          sn=sin(tmp)
          sums_omega_asymp(i_tau,1,1)=sums_omega_asymp(i_tau,1,1)
     &                                 +cs/w
          sums_omega_asymp(i_tau,1,2)=sums_omega_asymp(i_tau,1,2)
     &                                 +sn/w
          sums_omega_asymp(i_tau,2,1)=sums_omega_asymp(i_tau,2,1)
     &                                 +cs/w2
          sums_omega_asymp(i_tau,3,2)=sums_omega_asymp(i_tau,3,2)
     &                                 +sn/w3
        enddo
        sums_omega_asymp(i_tau,1,2)=betta_t/4.d0
     &                               -sums_omega_asymp(i_tau,1,2)
        if(i_tau==0) then
          tmp=log(4.d0)
        else
          tmp=log((1.d0-cos(2.d0*fr*t))/(1.d0-cos(fr*t)))
        endif
        sums_omega_asymp(i_tau,1,1)=betta_t/2.d0/pi*tmp
     &                               -sums_omega_asymp(i_tau,1,1)
        sums_omega_asymp(i_tau,2,1)=betta_t**2/4.d0*(0.5d0-t/betta_t)
     &                               -sums_omega_asymp(i_tau,2,1)
        sums_omega_asymp(i_tau,3,2)=betta_t**2*t/8.d0*(1.d0-t/betta_t)
     &                               -sums_omega_asymp(i_tau,3,2)
        do i=n_omega_maz+1,n_omega_maz+1000000
          w=(i+i+1.d0)*fr
          w2=w*w
          w3=w2*w
          w4=w2*w2
          w5=w3*w2
          w6=w3*w3
          tmp=w*t
          cs=cos(tmp)
          sn=sin(tmp)
          sums_omega_asymp(i_tau,3,1)=sums_omega_asymp(i_tau,3,1)
     &                                 +cs/w3
          sums_omega_asymp(i_tau,4,1)=sums_omega_asymp(i_tau,4,1)
     &                                 +cs/w4
          sums_omega_asymp(i_tau,5,1)=sums_omega_asymp(i_tau,5,1)
     &                                 +cs/w5
          sums_omega_asymp(i_tau,6,1)=sums_omega_asymp(i_tau,6,1)
     &                                 +cs/w6
          sums_omega_asymp(i_tau,2,2)=sums_omega_asymp(i_tau,2,2)
     &                                 +sn/w2
          sums_omega_asymp(i_tau,4,2)=sums_omega_asymp(i_tau,4,2)
     &                                 +sn/w4
          sums_omega_asymp(i_tau,5,2)=sums_omega_asymp(i_tau,5,2)
     &                                 +sn/w5
          sums_omega_asymp(i_tau,6,2)=sums_omega_asymp(i_tau,6,2)
     &                                 +sn/w6
          aa=max(abs(cs),abs(sn))/w2
          if(aa<1.d-20) exit
        enddo
      enddo
	if(nproc/=1) then
#ifdef MPI
	  call dgop(sums_omega_asymp,6*2*(n_tau/2+1),'  +',
     &	          MPI_COMM_WORLD)
#endif     
        sums_omega_asymp=(sums_omega_asymp*nproc_tau)/nproc
      endif
      sums_omega_asymp=sums_omega_asymp/betta_t
      end