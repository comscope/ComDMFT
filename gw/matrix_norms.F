      subroutine matrix_norms(a,n,n0,norm_1,norm_inf,ndiag_ratio)
c     norm_1 = Max Sum_i |A_ij| among columns
c     norm_inf = Max Sum_j |A_ij| among rows
c     ndiag_ratio = Sum(nondiag)|A_ij| / Sum(all)|A_ij|/(N-1)
	implicit none
	integer, intent(in) :: n,n0
	complex*16, intent(in) :: a(n0,n0)
	real*8, intent(out) :: norm_1,norm_inf,ndiag_ratio
	integer :: i,j
	real*8 :: d,dn,t,tn
	d=0.d0
	dn=0.d0
	norm_1=0.d0
	do j=1,n
	  t=0.d0
	  tn=0.d0
	  do i=1,n
	    t=t+abs(a(i,j))
	    if(i/=j) tn=tn+abs(a(i,j))
	  enddo
	  d=d+t
	  dn=dn+tn
	  norm_1=max(norm_1,t)
	enddo
	ndiag_ratio=0.d0
	if(n>1) ndiag_ratio=dn/d/(n-1)
	norm_inf=0.d0
	do i=1,n
	  t=0.d0
	  do j=1,n
	    t=t+abs(a(i,j))
	  enddo
	  norm_inf=max(norm_inf,t)
	enddo
      end