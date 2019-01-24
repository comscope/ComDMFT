      subroutine int_tau_sigc_gx(e,sig_tau,ex,n,n0)
c ------ Performs the integral of Sigc(-t)*Gx(t) -----------------------
	use manager_mod
	use parallel_mod
      implicit none
      integer, intent(in) :: n,n0
	real*8, intent(in) :: sig_tau(n0,n0,2,ndim3_tau),ex(n0)
	real*8, intent(out) :: e
	integer :: i,i_tau,ind_tau
	real*8 :: gt,gbt,green0,gx1,gx2
	real*8, allocatable :: f(:)
	allocate(f(ndim3_tau))
c ------------- AA nad BB parts ------------------------------------
      f=0.d0
      do ind_tau=1,ndim3_tau
        i_tau=me3_tau*ndim3_tau+ind_tau-1
	  do i=1,n
	    gt=green0(ex(i)-chem_pot,tau_mesh(i_tau))
	    gbt=green0(ex(i)-chem_pot,tau_mesh(n_tau-i_tau))
	    gx1=gt-gbt
	    gx2=gt+gbt
	    f(ind_tau)=f(ind_tau)+sig_tau(i,i,1,ind_tau)*gx1
     &	                       -sig_tau(i,i,2,ind_tau)*gx2
        enddo
      enddo
      e=dot_product(f,ab_nu_from_tau(:,0,1))/4.d0
	deallocate(f)
      end