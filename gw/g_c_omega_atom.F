      subroutine g_c_omega_atom(gx_omega,sigc,gc_omega)
      use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	complex*16, intent(in) :: gx_omega(maxntle,nrel*maxb+1),
     &	                        sigc(maxntle,maxntle,nrel*maxb+1)
	complex*16, intent(out) :: gc_omega(maxntle,maxntle,nrel*maxb+1)
	integer :: i,info,li,l,n,j
	integer, allocatable :: ipiv(:)
	complex*16, allocatable :: gs(:,:)
	allocate(gs(maxntle,maxntle))
	allocate(ipiv(maxntle))
	gc_omega=(0.d0,0.d0)
      do li=1,nrel*lmb(1)+1
	  if(irel.le.1) then
	    l=li-1
	  else
	    l=li/2
	  endif
        n=ntle(l,1)
c ---- GS = G_x * Sig_c ----------------------------------------------
        do j=1,n
          do i=1,n
            gs(i,j)=gx_omega(i,li)*sigc(i,j,li)
          enddo
        enddo
        do j=1,n
          do i=1,n
            gc_omega(i,j,li)=gs(i,j)*gx_omega(j,li)
          enddo
        enddo
        gs=-gs
	  do i=1,n
	    gs(i,i)=(1.d0,0.d0)+gs(i,i)
	  enddo
*     Compute the LU factorization of 1-GS.
*
        CALL ZGETRF(n,n,gs,maxntle,ipiv,INFO)
        IF( INFO.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
          CALL ZGETRS('n',n,n,gs,maxntle,ipiv,gc_omega(1,1,li),maxntle,
     &                info)
        END IF
      enddo
	deallocate(gs,ipiv)
      end