      subroutine pade_fermion_chain(ff,n,n0,ap,nc)
	use manager_mod
      implicit none
      integer, intent(in) :: n,n0,nc
	complex*16, intent(in) :: ff(n0,n0,0:n_omega)
	complex*16, intent(out) :: ap(0:nc,n0,n0)
	integer :: nn,i,j,ii,jj
	real*8 :: sm
	complex*16, allocatable :: f(:),g(:),zi(:)
	nn=nc
	allocate(f(0:nn))
	allocate(g(0:nn))
	allocate(zi(0:nn))
	zi=w_omega(0:nn)*(0.d0,1.d0)
	ap=(0.d0,0.d0)
	do jj=1,n
	  do ii=1,n
	    sm=sum(abs(ff(ii,jj,:)))
	    if(sm<1.d-12) cycle
	    ap(0,ii,jj)=ff(ii,jj,0)
	    f=ff(ii,jj,0:nn)
	    do j=1,nn
	      do i=j,nn
	        g(i)=(ap(j-1,ii,jj)-f(i))/(zi(i)-zi(j-1))/f(i)
	      enddo
	      ap(j,ii,jj)=g(j)
	      f=g
	    enddo
	  enddo
	enddo
	deallocate(f,g,zi)
      end
      
      
      subroutine ac_fermion_chain(ff0,ff0_cheb,n,n0,ap,nc)
	use manager_mod
	use units_mod
      implicit none
      integer, intent(in) :: n,n0,nc
	real*8, intent(in) :: ff0(n,n,2,0:n_omega),
     &	                    ff0_cheb(n,n,2,0:n_omega_asy-2)
	complex*16, intent(out) :: ap(0:nc,n0,n0)
	integer :: nn,i,j,ii,jj
	real*8 :: sm,pi2b
	real*8, allocatable :: w_out(:)
	complex*16, allocatable :: f(:),g(:),zi(:),ff(:,:,:)
	pi2b=(pi+pi)/betta_t
	nn=nc
	allocate(ff(n,n,0:nn))
	allocate(w_out(0:nn))
	do i=0,min(nn,n_omega_exa)
	  w_out(i)=w_omega(i)
	enddo
	do i=min(nn,n_omega_exa)+1,nn
	  w_out(i)=w_out(i-1)+pi2b
	enddo
	call interp_g_omega_c(ff0,ff0_cheb,ff,w_out,nn,n,0)
	allocate(f(0:nn))
	allocate(g(0:nn))
	allocate(zi(0:nn))
	zi=w_out(0:nn)*(0.d0,1.d0)
	ap=(0.d0,0.d0)
	do jj=1,n
	  do ii=1,n
	    sm=sum(abs(ff(ii,jj,:)))
	    if(sm<1.d-12) cycle
	    ap(0,ii,jj)=ff(ii,jj,0)
	    f=ff(ii,jj,0:nn)
	    do j=1,nn
	      do i=j,nn
	        g(i)=(ap(j-1,ii,jj)-f(i))/(zi(i)-zi(j-1))/f(i)
	      enddo
	      ap(j,ii,jj)=g(j)
	      f=g
	    enddo
	  enddo
	enddo
	deallocate(f,g,zi,ff,w_out)
      end