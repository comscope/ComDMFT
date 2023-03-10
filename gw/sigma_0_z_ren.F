! as input, tmp is the i_omega=0,1 elements of sig_c_omega
! as output, contains the omega=0 value of the self-energy, and the omega-derivative of the self-energy
! The tmp(:,:,1) elements contain sigma(omega=0)
! The tmp(:,:,2) elements contain the derivative dsigma/domega
! tmp0 is the diagonal elements of sig_c_omega, and is unchanged
! n is n_bnd(k,ispin)
! n0 is nbndf, the maximum value of n_bnd
! The only matrix elements that are touched are between (1:n), not (1:n0)
      subroutine sigma_0_z_ren(n,n0,tmp,tmp0)

	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod

      implicit none
	integer, intent(in) :: n,n0
	complex*16, intent(in) :: tmp0(n0,0:n_omega)
	complex*16, intent(inout) :: tmp(n0,n0,2)

	integer :: i,j,iom,kom
	real*8 :: hh
	complex*16 :: s1,s2,t1,t2
        kom=-1


	if(z_fit==1) then ! z_fit is hard-coded to 1

	  do j=1,n

! kom is the omega-index of first value of tmp0 satisfying (imag(tmp0(j,iom))<0.d0)
	    do iom=0,n_omega
	      if(imag(tmp0(j,iom))<0.d0) then
	        kom=iom
	        goto 1
	      endif
	    enddo
1         continue

          hh=2.d0*pi/betta_t*(2*kom+1)

! special work on the diagonal
	    tmp(j,j,1)=0.5d0*(tmp0(j,kom)+conjg(tmp0(j,kom)))
	    tmp(j,j,2)=(1.d0,0.d0)
     &	            +(0.d0,1.d0)*(tmp0(j,kom)-conjg(tmp0(j,kom)))/hh

	    do i=j+1,n


! The (:,:,2) elements contain the derivative dsigma/domega
	      if(kom==0) then
	        tmp(i,j,2)=(0.d0,1.d0)*(tmp(i,j,1)-conjg(tmp(j,i,1)))/hh
	      else if(kom>0) then
	        tmp(i,j,2)=(0.d0,0.d0)
	      endif

! The (:,:,1) elements contain sigma(omega=0)
	      tmp(i,j,1)=0.5d0*(tmp(i,j,1)+conjg(tmp(j,i,1)))

! Calculate the other half of the matrices.
	      tmp(j,i,1)=conjg(tmp(i,j,1))
	      tmp(j,i,2)=conjg(tmp(i,j,2))

	    enddo ! i=j+1,n

	  enddo

	else if(z_fit==2) then
	  do j=1,n
	    do i=1,j
	      s1=0.5d0*(tmp(i,j,1)+conjg(tmp(j,i,1)))
	      t1=0.5d0*(tmp(i,j,1)-conjg(tmp(j,i,1)))
	      s2=0.5d0*(tmp(i,j,2)+conjg(tmp(j,i,2)))
	      t2=0.5d0*(tmp(i,j,2)-conjg(tmp(j,i,2)))
	      tmp(i,j,1)=(9.d0*s1-s2)/8.d0
	      tmp(i,j,2)=(0.d0,1.d0)*betta_t*(27.d0*t1-t2)/24.d0/pi
	      if(i==j) then
	        tmp(j,i,2)=(1.d0,0.d0)+tmp(i,j,2)
	      else if(i/=j) then
	        tmp(j,i,1)=conjg(tmp(i,j,1))
	        tmp(j,i,2)=conjg(tmp(i,j,2))
	      endif
	    enddo
	  enddo
	endif
      end
