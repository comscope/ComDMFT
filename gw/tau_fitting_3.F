      subroutine tau_fitting_3(e0,a0,gt,nf_fit,err,w0,w1,n)
	use manager_mod
      implicit none
      integer, intent(in) :: nf_fit,n
	real*8, intent(in) :: gt(0:n_tau),w0(3),w1(3)
	real*8, intent(out) :: e0(nf_fit),a0(nf_fit),err
	integer :: i_tau,i1,i2,i3
	real*8 :: d,s,green0,hh1,hh2,hh3,g,a,gg,c1,c2,c3,e1,e2,e3,b,a1,b1,
     &	        st,tt,t
	real*8, allocatable :: g1(:),g2(:),g3(:,:)
	allocate(g1(0:n_tau))
	allocate(g2(0:n_tau))
	allocate(g3(0:n_tau,0:n))
	hh1=(w1(1)-w0(1))/n
	hh2=(w1(2)-w0(2))/n
	hh3=(w1(3)-w0(3))/n
	    do i3=0,n
	      e3=w0(3)+hh3*i3
	      do i_tau=0,n_tau
	        g3(i_tau,i3)=green0(e3,tau_mesh(i_tau))
	      enddo
	    enddo
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
	    do i3=0,n
	      e3=w0(3)+hh3*i3
	      if(e3-e2<2*hh3) cycle
c	      do i_tau=0,n_tau
c	        g3(i_tau)=green0(e3,tau_mesh(i_tau))
c	      enddo
	      gg=g2(0)-g3(0,i3)
	      a=(gt(0)-g3(0,i3))/gg
	      b=(g1(0)-g3(0,i3))/gg
	      a1=1.d0-a
	      b1=1.d0-b
	      st=0.d0
	      tt=0.d0
	      do i_tau=1,n_tau-1
	        s=gt(i_tau)-a*g2(i_tau)-a1*g3(i_tau,i3)
	        t=-g1(i_tau)+b*g2(i_tau)+b1*g3(i_tau,i3)
	        st=st+s*t
	        tt=tt+t*t
	      enddo
	      c1=-st/tt
	      c2=a-b*c1
	      c3=1.d0-c2-c1
	      s=0.d0
	      do i_tau=0,n_tau
	        g=c1*g1(i_tau)+c2*g2(i_tau)+c3*g3(i_tau,i3)
	        s=s+(gt(i_tau)-g)**2
	      enddo
	      if(s<d) then
	        d=s
	        e0(1)=e1
	        e0(2)=e2
	        e0(3)=e3
	        a0(1)=c1
	        a0(2)=c2
	        a0(3)=c3
	      endif
	    enddo
	  enddo
	enddo
	deallocate(g1,g2,g3)
	err=sqrt(d/(n_tau+1))
      end