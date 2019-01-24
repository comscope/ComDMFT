      double precision function green0_boson_tau(kap,tau)
c     B(t)=1/2/kap * (e^(-kap*t)+e^(-kap*(b-t)))/(1-e^(-kap*b))      
	use manager_mod
      implicit none
	real*8, intent(in) :: kap,tau
	real*8 :: a,b,t1,t2,t3
	t1=kap*tau
	t2=kap*betta_t
	t3=kap*(betta_t-tau)
	a=0.d0
	if(t1<100.d0) a=exp(-t1)
	if(t3<100.d0) a=a+exp(-t3)
	b=1.d0
	if(t2<100.d0) b=b-exp(-t2)
	green0_boson_tau=0.5d0*a/b/kap
      end