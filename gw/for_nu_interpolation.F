      subroutine for_nu_interpolation
      use manager_mod
      use units_mod
      implicit none
	real*8 :: w,w1
	allocate(asy_nu_interp(2,2,2))
      w=w_nu(n_nu)
      w1=w_nu(n_nu-1)
      asy_nu_interp(1,1,1)=w1**4/(w1-w)/(w1+w)
      asy_nu_interp(2,1,1)=w**4/(w-w1)/(w1+w)
      asy_nu_interp(1,2,1)=-w1**4*w**2/(w1-w)/(w1+w)
      asy_nu_interp(2,2,1)=w1**2*w**4/(w1-w)/(w1+w)
      asy_nu_interp(1,1,2)=w1**5/(w1-w)/(w1+w)
      asy_nu_interp(2,1,2)=w**5/(w-w1)/(w1+w)
      asy_nu_interp(1,2,2)=-w1**5*w**2/(w1-w)/(w1+w)
      asy_nu_interp(2,2,2)=w1**2*w**5/(w1-w)/(w1+w)
      end