      subroutine g_qp_tau(ispin,ind_k,gqp,ta,n,chem)
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: ispin,ind_k,n
      real*8, intent(in) :: chem,ta
      complex*16, intent(out) :: gqp(n,n)
	integer :: k,i,ib,jb
	real*8 :: de,green0
	real*8 :: g_diag
	k=n3_mpi_k(me3_k+1)+ind_k
      gqp=(0.d0,0.d0)
c ----- Full QP Green's function (OMEGA) -----------------------------
     	do i=1,n
        de=e_qp(i,k,ispin)-chem
	  g_diag=green0(de,ta)
	  do jb=1,n
	    do ib=1,n
	      gqp(ib,jb)=gqp(ib,jb)+q_qp(ib,i,ind_k,ispin)*g_diag
     &                           *conjg(q_qp(jb,i,ind_k,ispin))
	    enddo
        enddo
      enddo
      end
