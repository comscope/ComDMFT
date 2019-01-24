      subroutine bos_sums
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
	integer :: ind_tau,i_nu,i,i_tau
	real*8 :: fr,c1,cs,sn,tmp,w,t,w2,w3,w4,w5,w6,aa,w8,pi_bet,h1,x4,
     &	        x6,x8,nu1,nu2,nu3,dt,a1,a2,a3,b1,b2,b3,c2,c3,v,v2,v4,v6,
     &          v8
	pi_bet=pi/betta_t
	h1=pi_bet+pi_bet
c
c ------- Here sums for asymptotic contribution to TAU-functions ---
c
      sums_nu_asymp=0.d0
	do ind_tau=1,ndim_tau(me+1)
	  i_tau=n_mpi_tau(me+1)+ind_tau-1
        t=tau_mesh(i_tau)
        fr=2.d0*pi/betta_t
c ---- We use some known sums here for cos/w; cos/w2; ...----       
        do i=1,n_nu_maz
          w=i*fr
          w2=w*w
          tmp=w*t
          cs=cos(tmp)
          sums_nu_asymp(i_tau,1)=sums_nu_asymp(i_tau,1)+cs/w
          sums_nu_asymp(i_tau,2)=sums_nu_asymp(i_tau,2)+cs/w2
        enddo
        if(i_tau==0) then
          tmp=log(1.d0/(1.d0-cos(fr*1.d-10)))/2.d0 !it is fro tau=1.d-10
        else
          tmp=log(1.d0/(1.d0-cos(fr*t)))/2.d0
        endif
        sums_nu_asymp(i_tau,1)=betta_t/2.d0/pi*tmp
     &                           -sums_nu_asymp(i_tau,1)
        sums_nu_asymp(i_tau,2)=0.25d0*(t*t-betta_t*t+betta_t**2/6.d0)
     &                           -sums_nu_asymp(i_tau,2)
        do i=n_nu_maz+1,n_nu_maz+1000000
          w=i*fr
          w2=w*w
          w3=w2*w
          w4=w2*w2
          w5=w3*w2
          w6=w3*w3
          w8=w4*w4
          tmp=w*t
          cs=cos(tmp)
          sn=sin(tmp)
          sums_nu_asymp(i_tau,3)=sums_nu_asymp(i_tau,3)+cs/w3
          sums_nu_asymp(i_tau,4)=sums_nu_asymp(i_tau,4)+cs/w4
          sums_nu_asymp(i_tau,5)=sums_nu_asymp(i_tau,5)+cs/w5
          sums_nu_asymp(i_tau,6)=sums_nu_asymp(i_tau,6)+cs/w6
          sums_nu_asymp(i_tau,8)=sums_nu_asymp(i_tau,8)+cs/w8
          aa=abs(cs)/w3
          if(aa<1.d-20) exit
        enddo
      enddo 
	if(nproc_tau/=1) then
#ifdef MPI
	  call dgop(sums_nu_asymp,8*(n_tau/2+1),'  +',MPI_COMM_WORLD)
#endif	  
	endif
      sums_nu_asymp=sums_nu_asymp/betta_t
c ---- Here are Nu-Sums for Internal Energy evaluation ---------
      allocate(et_sum_nu(0:n_nu))
	et_sum_nu=0.d0	    
c ---------------------- Exact Part ------------------------------
	et_sum_nu(0)=0.5d0	  
	do i_nu=1,n_nu_exa
	  et_sum_nu(i_nu)=1.d0
	enddo
c ----------------- Asymptotic Part ------------------------------	
	do i=n_nu_exa+1,n_nu-1,2
	  nu1=w_nu(i-1)
	  nu2=w_nu(i)
	  nu3=w_nu(i+1)
c ---------------------- COS transformation -----------------------	    
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
	  do i_nu=mats_num_nu(i-1)+1,mats_num_nu(i+1)
	    v=h1*i_nu
	    v2=v*v
	    v4=v2*v2
	    v6=v4*v2
	    v8=v4*v4
	    x4=x4+1.d0/v4
	    x6=x6+1.d0/v6
	    x8=x8+1.d0/v8
	  enddo
	  et_sum_nu(i-1)=et_sum_nu(i-1)+x4*a1+x6*b1+x8*c1
	  et_sum_nu(i)=et_sum_nu(i)+x4*a2+x6*b2+x8*c2
	  et_sum_nu(i+1)=et_sum_nu(i+1)+x4*a3+x6*b3+x8*c3
	enddo
	et_sum_nu=et_sum_nu/betta_t	
      end