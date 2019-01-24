      subroutine ferm_interp_omega(n,ab_in,ab_cheb,omega,ab_out,key)
c ---- Interpolation A,B (fermionic) onto Matsubara's OMEGA --
c --- Input and output arrays are in compact R*8 form -----------------
c     key = 0 - for Gx+Gc
c     key = 1 - only for Gc
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n,key
	real*8, intent(in) :: ab_in(n,2,0:n_omega),omega,
     &                      ab_cheb(n,2,0:n_omega_asy-2)
	real*8, intent(out) :: ab_out(n,2)
	integer :: i,oma_mats,i1,j
	real*8 :: x,x_cheb_from_omega,v2,v4,pi8,bt,v,pi2b,v3,v5,v6
	real*8, allocatable :: cheb(:)
	pi8=8.d0*pi
	bt=pi/betta_t
	pi2b=bt+bt
	v=abs(omega)
      if(v<=w_omega(n_omega_exa)+1.d-12) then
        oma_mats=idnint((v-bt)/pi2b)
	  v2=bt+oma_mats*pi2b
	  if(abs(v-v2)<1.d-13) then
          ab_out=ab_in(:,:,oma_mats)
        else
	    i1=oma_mats-2
	    i1=max(i1,0)
	    i1=min(i1,n_omega_exa-5)
	    do j=1,2
	      do i=1,n
	        call fint(ab_in(i,j,i1:i1+5),w_omega_d(i1),6,v,v4)
	        ab_out(i,j)=v4
	      enddo
	    enddo
        endif
      else if(abs(v-w_omega(n_omega))<1.d-12) then
        ab_out=ab_in(:,:,n_omega)
      else if(v<w_omega(n_omega)) then
	  allocate(cheb(0:n_omega_asy-2))
        x=x_cheb_from_omega(v)
	  call dcheb_recurr(x,cheb,n_omega_asy-2)
	  ab_out=0.d0
	  do i=0,n_omega_asy-2
	    ab_out=ab_out+ab_cheb(:,:,i)*cheb(i)
	  enddo
	  deallocate(cheb)
      else
        v2=v*v
        v3=v2*v
        v4=v2*v2
        v5=v4*v
        v6=v4*v2
        if(key==0) then
	    ab_out(:,1)=
     &	    (asy_omega_interp(1,1,1)*ab_in(:,1,n_omega-1)
     &	    +asy_omega_interp(2,1,1)*ab_in(:,1,n_omega))/v2
     &      +(asy_omega_interp(1,2,1)*ab_in(:,1,n_omega-1)
     &	    +asy_omega_interp(2,2,1)*ab_in(:,1,n_omega))/v4
	    ab_out(:,2)=
     &	    (asy_omega_interp(1,1,2)*ab_in(:,2,n_omega-1)
     &      +asy_omega_interp(2,1,2)*ab_in(:,2,n_omega))/v
     &      +(asy_omega_interp(1,2,2)*ab_in(:,2,n_omega-1)
     &      +asy_omega_interp(2,2,2)*ab_in(:,2,n_omega))/v3
        else if(key==1) then
	    ab_out(:,1)=
     &	    (asy_omega_interp_corr(1,1,1)*ab_in(:,1,n_omega-1)
     &	    +asy_omega_interp_corr(2,1,1)*ab_in(:,1,n_omega))/v4
     &      +(asy_omega_interp_corr(1,2,1)*ab_in(:,1,n_omega-1)
     &	    +asy_omega_interp_corr(2,2,1)*ab_in(:,1,n_omega))/v6
	    ab_out(:,2)=
     &	    (asy_omega_interp_corr(1,1,2)*ab_in(:,2,n_omega-1)
     &      +asy_omega_interp_corr(2,1,2)*ab_in(:,2,n_omega))/v3
     &      +(asy_omega_interp_corr(1,2,2)*ab_in(:,2,n_omega-1)
     &      +asy_omega_interp_corr(2,2,2)*ab_in(:,2,n_omega))/v5
        endif
	endif
      end