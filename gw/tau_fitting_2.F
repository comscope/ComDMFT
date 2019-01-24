      subroutine tau_fitting_2(e0,a0,gt,nf_fit,err,w0,w1,n)
	use manager_mod
      implicit none
      integer, intent(in) :: nf_fit,n
	real*8, intent(in) :: gt(0:n_tau),w0(2),w1(2)
	real*8, intent(out) :: e0(nf_fit),a0(nf_fit),err
	integer :: i_tau,i1,i2
	real*8 :: d,s,green0,hh1,hh2,g,c1,c2,e1,e2,g0,gb
	real*8, allocatable :: g1(:),g2(:)
	allocate(g1(0:n_tau))
	allocate(g2(0:n_tau))
	g0=gt(0)
	gb=gt(n_tau)
	hh1=(w1(1)-w0(1))/n
	hh2=(w1(2)-w0(2))/n
	d=1.d5
	do i1=0,n
	  e1=w0(1)+hh1*i1
	  do i_tau=0,n_tau
	    g1(i_tau)=green0(e1,tau_mesh(i_tau))
	  enddo
	  do i2=0,n
	    e2=w0(2)+hh2*i2
	    if(e2-e1<2*hh2) cycle
	    do i_tau=0,n_tau
	      g2(i_tau)=green0(e2,tau_mesh(i_tau))
	    enddo
	    if(abs(e1*betta_t)<1.d-10.and.abs(e2*betta_t)<1.d-10) then
	      c1=1.d0
	    else if(abs(g0)>abs(gb)) then
	      c1=(g0-g2(0))/(g1(0)-g2(0))
	    else
	      c1=(gb-g2(n_tau))/(g1(n_tau)-g2(n_tau))
	    endif
	    c2=1.d0-c1
	    s=0.d0
	    do i_tau=0,n_tau
	      g=c1*g1(i_tau)+c2*g2(i_tau)
	      s=s+(gt(i_tau)-g)**2
	    enddo
	    if(s<d) then
	      d=s
	      e0(1)=e1
	      e0(2)=e2
	      a0(1)=c1
	      a0(2)=c2
	    endif
	  enddo
	enddo
	deallocate(g1,g2)
	err=sqrt(d/(n_tau+1))
      end