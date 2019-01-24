      subroutine sum_maz_sig_gx(e,sigc_omega,ex,n,n1)
c - Performs the summation of Re{Sigma(w)Gc(w)} over Matsubara frequences
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n1
	real*8, intent(in) :: sigc_omega(n1,n1,2,ndim3_omega),ex(n)
	real*8, intent(out) :: e
	integer :: i_omega,i,ind_omega
	real*8 :: x,y,w,et,g1,g2
      e=0.d0
	do ind_omega=1,ndim3_omega
	  i_omega=me3_tau*ndim3_omega+ind_omega-1
	  w=w_omega(i_omega)
	  et=0.d0
	  do i=1,n
	    x=chem_pot-ex(i)
	    y=w*w+x*x
	    g1=x/y
	    g2=-w/y
          et=et+sigc_omega(i,i,1,ind_omega)*g1
     &         -sigc_omega(i,i,2,ind_omega)*g2
        enddo
        e=e+et_sum_sc_gx(i_omega)*et
      enddo
      end