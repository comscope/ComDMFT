      subroutine asympt_sums_omega
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
	integer :: ind_tau,i_tau,k,n
	real*8 :: su(3,2),ber(2),ber0(2),cs,sn,t,x,x0,pi_bet,a,a0,pi2,xp,
     &	        tmp,pi0,k2,k2p
      pi0=3.14159265358979323846d0
	pi_bet=pi0/betta_t
	a=pi_bet+pi_bet
	a0=pi_bet
	pi2=pi0+pi0
c
c ------- Here sums for asymptotic contribution to TAU-functions ---
c
      allocate(sum_omega_asy(0:n_tau/2,3,2))
      sum_omega_asy=0.d0
	do ind_tau=1,ndim3_tau
	  i_tau=me3_tau*ndim3_tau+ind_tau-1
        t=tau_mesh(i_tau)
        x=a*t
        x0=a0*t
c ---- Get Bernoulli polynomials ---
        xp=x/pi2
        ber(1)=xp-0.5d0
        ber(2)=xp**2-xp+1.d0/6.d0
        xp=x0/pi2
        ber0(1)=xp-0.5d0
        ber0(2)=xp**2-xp+1.d0/6.d0
c ---- Get partial sums ---   
        su=0.d0    
        do k=0,n_omega_maz
          k2p=k+k+1.d0
          k2=k2p*k2p
          tmp=k2p*x0
          cs=cos(tmp)/k2
          sn=sin(tmp)/k2p
          su(1,1)=su(1,1)+cs
          su(1,2)=su(1,2)+sn
        enddo   
        do k=n_omega_maz+1,n_omega_maz*128
          k2p=k+k+1.d0
          k2=k2p*k2p
          tmp=k2p*x0
          cs=cos(tmp)/k2
          sn=sin(tmp)/k2p
          cs=cs/k2
          sn=sn/k2
          su(2,1)=su(2,1)+cs
          su(2,2)=su(2,2)+sn
          cs=cs/k2
          sn=sn/k2
          su(3,1)=su(3,1)+cs
          su(3,2)=su(3,2)+sn
        enddo
        cs=pi0**2
        sn=-pi0
        sum_omega_asy(i_tau,1,1)=
     &    (cs*(ber0(2)-ber(2)/4.d0)-su(1,1))/a0**2
        sum_omega_asy(i_tau,1,2)=
     &    (sn*(ber0(1)-ber(1)/2.d0)-su(1,2))/a0
        do n=2,3
          sum_omega_asy(i_tau,n,1)=su(n,1)/a0**(n+n)
          sum_omega_asy(i_tau,n,2)=su(n,2)/a0**(n+n-1)
        enddo
      enddo   !! over ind_tau
	if(nproc_tau/=1) call dgop(sum_omega_asy,6*(n_tau/2+1),'  +',
     &	                         comm_pnk)
      end