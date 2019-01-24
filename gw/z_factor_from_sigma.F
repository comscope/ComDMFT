      subroutine z_factor_from_sigma(sig_omega,n,n0,zfac)
c	This is for local Self Energies - we store them as complex matrices
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0
	complex*16, intent(in) :: sig_omega(n0,n0,ndim3_omega)
	complex*16, intent(out) :: zfac(n0,n0)
	integer :: lwork,info,i,j
	complex*16 :: iw
	integer, allocatable :: ipiv(:)
	complex*16, allocatable :: work(:)
c ------ Slopes of Sigma at Omega=0 -----------------------
      if(me3_tau==0) then
	  iw=dcmplx(0.d0,betta_t/2.d0/pi)
	  do j=1,n
	    do i=1,n
	      zfac(i,j)=iw*(sig_omega(i,j,1)-conjg(sig_omega(j,i,1))) 
		  if(i==j) zfac(i,j)=1.d0+zfac(i,j)
	    enddo
	  enddo
	endif
	if(nproc_tau/=1) then
	  call brdcst(msgdbl,zfac,16*n0*n0,0,comm_pnk)
	endif
c  Matrix ZFAC is inverse Z factor now ---------------------------
	lwork=3*n
	allocate(ipiv(n),work(lwork))
	call zgetrf(n,n,zfac,n0,ipiv,info)
	call zgetri(n,zfac,n0,ipiv,work,lwork,info)
	deallocate(ipiv,work)
      end