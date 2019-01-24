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
	if(z_fit==1) then
	  do j=1,n
	    do iom=0,n_omega
	      if(imag(tmp0(j,iom))<0.d0) then
	        kom=iom
	        goto 1
	      endif
	    enddo
1         continue
          hh=2.d0*pi/betta_t*(2*kom+1)
	    tmp(j,j,1)=0.5d0*(tmp0(j,kom)+conjg(tmp0(j,kom)))
	    tmp(j,j,2)=(1.d0,0.d0)
     &	            +(0.d0,1.d0)*(tmp0(j,kom)-conjg(tmp0(j,kom)))/hh
	    do i=j+1,n
	      if(kom==0) then
	        tmp(i,j,2)=(0.d0,1.d0)*(tmp(i,j,1)-conjg(tmp(j,i,1)))/hh
	      else if(kom>0) then
	        tmp(i,j,2)=(0.d0,0.d0)
	      endif
	      tmp(i,j,1)=0.5d0*(tmp(i,j,1)+conjg(tmp(j,i,1)))
	      tmp(j,i,1)=conjg(tmp(i,j,1))
	      tmp(j,i,2)=conjg(tmp(i,j,2))
	    enddo
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
