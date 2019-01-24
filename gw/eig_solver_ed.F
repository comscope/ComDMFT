      subroutine eig_solver_ed(n,nbas_full,ind_bas,ind_full,nbas,e,psi,
     &                         emm1,uat)
*****************************************************************************
*     Diagonalizes the Hamiltonian matrix for Exact Diagonalization Routine *
*****************************************************************************
c	nel - number of electrons
	use manager_mod
	use models_mod
	use parallel_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,nbas_full,nbas,ind_bas(nbas_max),
     &	                     ind_full(0:nbas_full-1)
	real*8, intent(in) :: emm1(n,n,nspin),uat(n,n,n,n,nspin,nspin)
	real*8, intent(out) :: e(nbas_max),psi(nbas_max,nbas_max)
	integer :: info,k,lwork,liwork,il,iu,ne
	real*8 :: vl,vu
	integer, allocatable :: iwork(:),ifail(:)
	real*8, allocatable :: work(:),a(:,:)
	do k=1,nbas !! over all distrib-s of n el-s over N One-P-states
        call h_action_atom(n,psi(:,k),ind_bas(k),nbas,ind_full,emm1,uat)
	enddo
	lwork=10*nbas+2*nbas**2
	liwork=3+5*nbas
      allocate(work(lwork),iwork(liwork),a(nbas,nbas),ifail(nbas))
c ---------- Eigen pairs of Hamiltonian -------------------------------
	a=psi(1:nbas,1:nbas)
      call DSYEVX_1('V','A','U',nbas,a,nbas,vl,vu,il,iu,1.d-35,ne,e,psi,
     &	            nbas_max,work,lwork,iwork,IFAIL, INFO )
	deallocate(work,iwork,a,ifail)
      end