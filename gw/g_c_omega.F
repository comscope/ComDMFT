      subroutine g_c_omega(key,gx_omega,gc_omega,ndim,n0)
c ------ Correlation part of Green's Function is calculated here -------------
c	key = 0  - For correlated part of G
c	key = 1  - For total G
c ----- In the beginning we have G_c = Sigma_c ---------------------------
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: ndim,n0,key
	complex*16, intent(in) :: gx_omega(ndim)
	complex*16, intent(inout) :: gc_omega(n0,n0)
	integer :: ipiv(ndim),i,j,info
	complex*16, allocatable :: a(:,:)
	allocate(a(ndim,ndim))
c ---- A = G_x * Sig_c ---------------------------------------------------
	do j=1,ndim
	  do i=1,ndim
	    a(i,j)=gx_omega(i)*gc_omega(i,j)
	  enddo
	enddo
	if(key==0) then
c ---- B = G_x * Sig_c * G_x  - Then it will be G_c -------------------------
	  do j=1,ndim
	    do i=1,ndim
	      gc_omega(i,j)=a(i,j)*gx_omega(j)
	    enddo
	  enddo
	else if(key==1) then
c ---- B = G_x  - Then it will be G_full -------------------------
        gc_omega=(0.d0,0.d0)
	  do j=1,ndim
	    gc_omega(j,j)=gx_omega(j)
	  enddo
	endif
c ---- A = I - G_x * Sig_c --------------------------------------------------
      call zscal(ndim**2,(-1.d0,0.d0),a,1)
      do i=1,ndim
	  a(i,i)=a(i,i)+(1.d0,0.d0)
	enddo
*     Compute the LU factorization of A.
*
      CALL ZGETRF(ndim,ndim,a,ndim,ipiv,INFO)
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
        CALL ZGETRS('n',ndim,ndim,a,ndim,ipiv,gc_omega,n0,info)
      END IF
	deallocate(a)
      end