      subroutine for_omega_interpolation
      use manager_mod
      use units_mod
      implicit none
	real*8 :: w,w1
c -------- For full (Gx+Gc) Green's function --------------------	
	allocate(asy_omega_interp(2,2,2))
      w=w_omega(n_omega)
      w1=w_omega(n_omega-1)
      asy_omega_interp(1,1,1)=w1**4/(w1-w)/(w1+w)
      asy_omega_interp(2,1,1)=w**4/(w-w1)/(w1+w)
      asy_omega_interp(1,2,1)=-w1**4*w**2/(w1-w)/(w1+w)
      asy_omega_interp(2,2,1)=w1**2*w**4/(w1-w)/(w1+w)
      asy_omega_interp(1,1,2)=w1**3/(w1-w)/(w1+w)
      asy_omega_interp(2,1,2)=w**3/(w-w1)/(w1+w)
      asy_omega_interp(1,2,2)=-w1**3*w**2/(w1-w)/(w1+w)
      asy_omega_interp(2,2,2)=w1**2*w**3/(w1-w)/(w1+w)
c -------- For correlation part (Gc) of Green's function --------
	allocate(asy_omega_interp_corr(2,2,2))
      asy_omega_interp_corr(1,1,1)=w1**6/(w1-w)/(w1+w)
      asy_omega_interp_corr(2,1,1)=w**6/(w-w1)/(w1+w)
      asy_omega_interp_corr(1,2,1)=-w1**6*w**2/(w1-w)/(w1+w)
      asy_omega_interp_corr(2,2,1)=w1**2*w**6/(w1-w)/(w1+w)
      asy_omega_interp_corr(1,1,2)=w1**5/(w1-w)/(w1+w)
      asy_omega_interp_corr(2,1,2)=w**5/(w-w1)/(w1+w)
      asy_omega_interp_corr(1,2,2)=-w1**5*w**2/(w1-w)/(w1+w)
      asy_omega_interp_corr(2,2,2)=w1**2*w**5/(w1-w)/(w1+w)
      end