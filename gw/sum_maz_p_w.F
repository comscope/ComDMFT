      subroutine sum_maz_p_w(e,p_nu,ws_nu,n)
c -------- Performs the summation of Re{P(v)W(v)} over Mazubara frequences
c     Exact Form here is: -1/beta (1/2 P(0)W(0) + Sum_{v>0}  P(v)W(v))
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
      integer, intent(in) :: n
	real*8, intent(in) :: p_nu(ndim3_nu,n,n),ws_nu(ndim3_nu,n,n)
	real*8, intent(out) :: e
	integer :: i_nu,ind_nu,i,j
	real*8 :: d
	e=0.d0
	do ind_nu=1,ndim3_nu
	  i_nu=me3_tau*ndim3_nu+ind_nu-1
	  d=0.d0
	  do j=1,n
	    do i=1,n
	      d=d+p_nu(ind_nu,i,j)*ws_nu(ind_nu,j,i)
	    enddo
	  enddo
	  e=e-d*et_sum_nu(i_nu)
	enddo
	if(nproc_tau/=1) call dgop(e,1,'  +',comm_pnk)
      end