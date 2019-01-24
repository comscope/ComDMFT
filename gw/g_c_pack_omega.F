      subroutine g_c_pack_omega(g0,d,dg,n,n0,n1,n2)
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: n,n0,n1,n2
	complex*16, intent(in) :: g0(n0,n0),d(n1,n1)
	real*8, intent(out) :: dg(n2,n2,2)
	integer :: i,info
	integer, allocatable :: ipiv(:)
	complex*16, allocatable :: a(:,:),tmp(:,:)
	allocate(a(n,n),ipiv(n),tmp(n,n))
c ---- A = G_0 * D ---------------------------------------------------
	call zgemm('n','n',n,n,n,(1.d0,0.d0),g0,n0,d,n1,(0.d0,0.d0),a,n)
	call zgemm('n','n',n,n,n,(1.d0,0.d0),a,n,g0,n0,(0.d0,0.d0),tmp,n)
	a=-a
      do i=1,n
	  a(i,i)=a(i,i)+(1.d0,0.d0)
	enddo
      CALL ZGETRF(n,n,a,n,ipiv,INFO)
      CALL ZGETRS('n',n,n,a,n,ipiv,tmp,n,info)
      call ferm_pack_omega(tmp,dg,n,n,n2)
	deallocate(a,ipiv,tmp)
      end