      subroutine transf_tau_from_omega(t,nt,tran)
c     Creates transformation matrices for OMEGA ---> TAU
c       They correspond to the transforms :
c     A(t) = 4/b Sum_w>0 cos(wt) A(w)
c     B(t) = 4/b Sum_w>0 sin(wt) B(w)
      use manager_mod
      use units_mod
      implicit none
	integer, intent(in) :: nt
      real*8, intent(in) :: t(nt)
      real*8, intent(out) :: tran(0:n_omega,nt,2)  !! COS  +  SIN parts
	integer :: i_tau,i_omega,i
	real*8 :: bt,pi2_bet,vt,c1,c2,c3,v,omega1,omega2,omega3,dt,a1,a2,
     &	        a3,b1,b2,b3,x2,x4,x6,x3,x5,x7,cs,v2,v4,v6,sn,v3,v5,v7
	bt=2.d0/betta_t
	pi2_bet=bt*pi
	bt=bt+bt
	tran=0.d0
	do i_tau=1,nt
c ---------------------- Exact Part ------------------------------	  
	  do i_omega=0,n_omega_exa
	    vt=w_omega(i_omega)*t(i_tau)
	    tran(i_omega,i_tau,1)=bt*cos(vt)
	    tran(i_omega,i_tau,2)=bt*sin(vt)
	  enddo
c ----------------- Asymptotic Part ------------------------------	  
	  do i=n_omega_exa+1,n_omega-1,2
	    omega1=w_omega(i-1)
	    omega2=w_omega(i)
	    omega3=w_omega(i+1)
c ---------------------- COS transformation -----------------------	    
	    dt=((1.d0/omega3**2-1.d0/omega2**2)/(omega2*omega3)**2
     &	     +(1.d0/omega2**2-1.d0/omega1**2)/(omega1*omega2)**2
     &	     +(1.d0/omega1**2-1.d0/omega3**2)/(omega1*omega3)**2)
     &       /(omega1*omega2*omega3)**2
	    a1=(1.d0/omega3**2-1.d0/omega2**2)/(omega2*omega3)**4/dt
	    a2=(1.d0/omega1**2-1.d0/omega3**2)/(omega1*omega3)**4/dt
	    a3=(1.d0/omega2**2-1.d0/omega1**2)/(omega1*omega2)**4/dt
	    b1=(1.d0/omega2**4-1.d0/omega3**4)/(omega2*omega3)**2/dt
	    b2=(1.d0/omega3**4-1.d0/omega1**4)/(omega1*omega3)**2/dt
	    b3=(1.d0/omega1**4-1.d0/omega2**4)/(omega1*omega2)**2/dt
	    c1=(1.d0/omega3**2-1.d0/omega2**2)/(omega2*omega3)**2/dt
	    c2=(1.d0/omega1**2-1.d0/omega3**2)/(omega1*omega3)**2/dt
	    c3=(1.d0/omega2**2-1.d0/omega1**2)/(omega1*omega2)**2/dt
	    x2=0.d0
	    x4=0.d0
	    x6=0.d0
	    do i_omega=mats_num_omega(i-1)+1,mats_num_omega(i+1)
	      v=pi2_bet/2.d0+pi2_bet*i_omega
	      cs=cos(v*t(i_tau))
	      v2=v*v
	      v4=v2*v2
	      v6=v4*v2
	      x2=x2+cs/v2
	      x4=x4+cs/v4
	      x6=x6+cs/v6
	    enddo
	    tran(i-1,i_tau,1)=tran(i-1,i_tau,1)+(x2*a1+x4*b1+x6*c1)*bt
	    tran(i,i_tau,1)=tran(i,i_tau,1)+(x2*a2+x4*b2+x6*c2)*bt
	    tran(i+1,i_tau,1)=tran(i+1,i_tau,1)+(x2*a3+x4*b3+x6*c3)*bt
c ---------------------- SIN transformation -----------------------	    
	    dt=((1.d0/omega3**2-1.d0/omega2**2)/(omega2*omega3)**2
     &	     +(1.d0/omega2**2-1.d0/omega1**2)/(omega1*omega2)**2
     &	     +(1.d0/omega1**2-1.d0/omega3**2)/(omega1*omega3)**2)
     &       /(omega1*omega2*omega3)**3
	    a1=(1.d0/omega3**2-1.d0/omega2**2)/(omega2*omega3)**5/dt
	    a2=(1.d0/omega1**2-1.d0/omega3**2)/(omega1*omega3)**5/dt
	    a3=(1.d0/omega2**2-1.d0/omega1**2)/(omega1*omega2)**5/dt
	    b1=(1.d0/omega2**4-1.d0/omega3**4)/(omega2*omega3)**3/dt
	    b2=(1.d0/omega3**4-1.d0/omega1**4)/(omega1*omega3)**3/dt
	    b3=(1.d0/omega1**4-1.d0/omega2**4)/(omega1*omega2)**3/dt
	    c1=(1.d0/omega3**2-1.d0/omega2**2)/(omega2*omega3)**3/dt
	    c2=(1.d0/omega1**2-1.d0/omega3**2)/(omega1*omega3)**3/dt
	    c3=(1.d0/omega2**2-1.d0/omega1**2)/(omega1*omega2)**3/dt
	    x3=0.d0
	    x5=0.d0
	    x7=0.d0
	    do i_omega=mats_num_omega(i-1)+1,mats_num_omega(i+1)
	      v=pi2_bet/2.d0+pi2_bet*i_omega
	      sn=sin(v*t(i_tau))
	      v2=v*v
	      v3=v2*v
	      v5=v3*v2
	      v7=v5*v2
	      x3=x3+sn/v3
	      x5=x5+sn/v5
	      x7=x7+sn/v7
	    enddo
	    tran(i-1,i_tau,2)=tran(i-1,i_tau,2)+(x3*a1+x5*b1+x7*c1)*bt
	    tran(i,i_tau,2)=tran(i,i_tau,2)+(x3*a2+x5*b2+x7*c2)*bt
	    tran(i+1,i_tau,2)=tran(i+1,i_tau,2)+(x3*a3+x5*b3+x7*c3)*bt
	  enddo  !! over i
	enddo  !! over i_tau
      end