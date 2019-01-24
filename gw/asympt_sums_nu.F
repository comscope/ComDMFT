      subroutine asympt_sums_nu
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
	integer :: ind_tau,i_tau,k,n
	real*8 :: cs,sn,tmp,t,pi_bet,ber(2),a,pi2,x,xp,pi0,su(3,2),k2
      pi0=3.1415926535897932384d0
	pi_bet=pi0/betta_t
	a=pi_bet+pi_bet
	pi2=pi0+pi0
c
c ------- Here sums for asymptotic contribution to TAU-functions ---
c
      allocate(sum_nu_asy(0:n_tau/2,3,2))
      sum_nu_asy=0.d0
	do ind_tau=1,ndim3_tau
	  i_tau=me3_tau*ndim3_tau+ind_tau-1
        t=tau_mesh(i_tau)
        x=a*t
c ---- Get Bernoulli polynomials ---
        xp=x/pi2
        ber(1)=xp-0.5d0
        ber(2)=xp**2-xp+1.d0/6.d0
c ---- Get partial sums ---
        su=0.d0
        do k=1,n_nu_maz
          k2=float(k)**2
          tmp=k*x
          cs=cos(tmp)/k2
          sn=sin(tmp)/k
          su(1,1)=su(1,1)+cs
          su(1,2)=su(1,2)+sn
        enddo
        do k=n_nu_maz+1,n_nu_maz*128
          k2=float(k)**2
          tmp=k*x
          cs=cos(tmp)/k2
          sn=sin(tmp)/k
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
        sum_nu_asy(i_tau,1,1)=(cs*ber(2)-su(1,1))/a**2
        sum_nu_asy(i_tau,1,2)=(sn*ber(1)-su(1,2))/a
        do n=2,3
          sum_nu_asy(i_tau,n,1)=su(n,1)/a**(n+n)
          sum_nu_asy(i_tau,n,2)=su(n,2)/a**(n+n-1)
        enddo
      enddo   !! over ind_tau
	if(nproc_tau/=1) call dgop(sum_nu_asy,6*(n_tau/2+1),'  +',
     &	                         comm_pnk)
      end
