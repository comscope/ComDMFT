      subroutine tau_fitting(e0,a0,gt,nf_fit,err)
	use manager_mod
      implicit none
      integer, intent(in) :: nf_fit
	real*8, intent(in) :: gt(0:n_tau)
	real*8, intent(out) :: e0(nf_fit),a0(nf_fit),err
	integer :: n,ie,i_tau,je
	real*8 :: d,s,green0,hh,ee,g,g0,gb,g1,g2,a,ee1,w0,w1
	w0=-0.5d0
	w1=1.d0
	n=600
	hh=(w1-w0)/n
	d=1.d5
	g0=gt(0)
	gb=gt(n_tau)
	do je=0,n
	  ee1=w0+hh*je
	  do ie=je+1,n
	    ee=w0+hh*ie
	    if(abs(ee*betta_t)<1.d-10.and.abs(ee1*betta_t)<1.d-10) then
	      a=1.d0
	    else if(abs(g0)>abs(gb)) then
	      g1=green0(ee1,0.d0)
	      g2=green0(ee,0.d0)
	      a=(g0-g2)/(g1-g2)
	    else
	      g1=green0(ee1,betta_t)
	      g2=green0(ee,betta_t)
	      a=(gb-g2)/(g1-g2)
	    endif
	    s=0.d0
	    do i_tau=0,n_tau
	      g1=green0(ee1,tau_mesh(i_tau))
	      g2=green0(ee,tau_mesh(i_tau))
	      g=a*g1+(1.d0-a)*g2
	      s=max(s,abs(gt(i_tau)-g))
	    enddo
	    if(s<d) then
	      d=s
	      e0(1)=ee1
	      e0(2)=ee
	      a0(1)=a
	      a0(2)=1.d0-a
	    endif
	  enddo
	enddo
	err=d
      end