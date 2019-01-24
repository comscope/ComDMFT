      double precision function x_cheb_from_tau(t)
      use manager_mod
      implicit none
      real*8, intent(in) :: t
      x_cheb_from_tau=(2.d0*t/betta_t)**(1.d0/3.d0)-1.d0
      end
      
      double precision function x_cheb_from_omega(w)
      use manager_mod
      implicit none
      real*8, intent(in) :: w
      real*8 :: w0,wm,a,b,c
      w0=sqrt(w_omega(n_omega_exa))
      wm=sqrt(omega_max)
      a=2.d0*w0*wm
      b=w0+wm
      c=wm-w0
      x_cheb_from_omega=(b-a/sqrt(w))/c
      end
      
      double precision function x_cheb_from_nu(v)
      use manager_mod
      implicit none
      real*8, intent(in) :: v
      real*8 :: v0,vm,a,b,c
      v0=sqrt(w_nu(n_nu_exa))
      vm=sqrt(nu_max)
      a=2.d0*v0*vm
      b=v0+vm
      c=vm-v0
      x_cheb_from_nu=(b-a/sqrt(v))/c
      end