      subroutine sum_maz_sig_gc(e,sigc_omega,gc_omega,n,n0,n1)
c - Performs the summation of Re{Sigma(w)Gc(w)} over Matsubara freq-s -
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: gc_omega(n0,n0,2,ndim3_omega),
     &	                    sigc_omega(n1,n1,2,ndim3_omega)
	real*8, intent(out) :: e
	integer :: i_omega,j,i,ind_omega
	real*8 :: de
	e=0.d0
	do j=1,n
	  do i=1,n
c --------- Exact contribution -------------------------------------
          de=0.d0
	    do ind_omega=1,ndim3_omega
	      i_omega=me3_tau*ndim3_omega+ind_omega-1
            de=de+et_sum_sc_gc(i_omega)
     &        *(sigc_omega(i,j,1,ind_omega)*gc_omega(i,j,1,ind_omega)
     &         -sigc_omega(i,j,2,ind_omega)*gc_omega(i,j,2,ind_omega))
          enddo
          if(i==j) e=e+de
          if(i/=j) e=e+de+de
        enddo
      enddo
      end