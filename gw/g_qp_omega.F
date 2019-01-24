      subroutine g_qp_omega(ispin,ind_k,gqp,w_n,n,chem)
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: ispin,ind_k,n
      real*8, intent(in) :: chem,w_n
      complex*16, intent(out) :: gqp(n,n)
	integer :: k,i,ib,jb
	real*8 :: de
	complex*16 :: g_diag
	k=n3_mpi_k(me3_k+1)+ind_k
      gqp=(0.d0,0.d0)
c ----- Full QP Green's function (OMEGA) -----------------------------
     	do i=1,n
        de=e_qp(i,k,ispin)-chem
	  g_diag=(1.d0,0.d0)/(dcmplx(0.d0,w_n)-de)
	  do jb=1,n
	    do ib=1,n
	      gqp(ib,jb)=gqp(ib,jb)+q_qp(ib,i,ind_k,ispin)*g_diag
     &                           *conjg(q_qp(jb,i,ind_k,ispin))
	    enddo
        enddo
      enddo
      end
