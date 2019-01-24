      subroutine eig_val_gener_solver(eps,n,n0,a,b,new,da,db)
      implicit none
	integer, intent(in) :: n,n0
	real*8, intent(in) :: eps
	complex*16, intent(inout) :: a(n0,n0),b(n0,n0)
	integer, intent(out) :: new
	real*8, intent(out) :: da(n),db(n)
      integer :: i,j
      real*8 :: s
	complex*16, allocatable :: z(:,:),a1(:,:)
      allocate(z(n,n),a1(n,n)) 
c  --------- Diagonalize Overlap Matrix ----------------------------
      do j=1,n
        do i=1,n
          z(i,j)=b(i,j)
        enddo
      enddo
      call eig_val_solver(n,n,z,db)
      new=0
      do i=1,n
        if(db(i)>eps) then
          new=new+1
          s=sqrt(db(i))
          b(:,new)=z(:,i)/s
        endif
      enddo
c ---------- Transform Hamiltonian to the optimized basis -----------
	call zgemm('c','n',new,n,n,(1.d0,0.d0),b,n,a,n,(0.d0,0.d0),a1,n) 
	call zgemm('n','n',new,new,n,(1.d0,0.d0),a1,n,b,n,(0.d0,0.d0),z,n)   
c  --------- Diagonalize Hamiltonian Matrix ----------------------------
      call eig_val_solver(new,n,z,da)
	call zgemm('n','n',n,new,new,(1.d0,0.d0),b,n,z,n,(0.d0,0.d0),a,n)
	deallocate(z,a1)
      end