      subroutine pade_fermion_chain_out(n,n0,ap,om,ff,nc)
	use manager_mod
	use units_mod
      implicit none
      integer, intent(in) :: n,n0,nc
      complex*16, intent(in) :: ap(0:nc,n0,n0),om
	complex*16, intent(out) :: ff(n0,n0)
	integer :: nn,j,ii,jj,i
	real*8 :: pi2b
	complex*16 :: c1
	real*8, allocatable :: w_out(:)
	complex*16, allocatable :: zi(:)
	pi2b=(pi+pi)/betta_t
	nn=nc
	allocate(zi(0:nn))
	allocate(w_out(0:nn))
	do i=0,min(nn,n_omega_exa)
	  w_out(i)=w_omega(i)
	enddo
	do i=min(nn,n_omega_exa)+1,nn
	  w_out(i)=w_out(i-1)+pi2b
	enddo
	zi=w_out(0:nn)*(0.d0,1.d0)
c --------------------------------------------------------------
	do jj=1,n
	  do ii=1,n
          c1=(1.d0,0.d0)
          do j=nn,1,-1
            c1=(1.d0,0.d0)+(om-zi(j-1))*ap(j,ii,jj)/c1
          enddo
          ff(ii,jj)=ap(0,ii,jj)/c1
        enddo
      enddo
	deallocate(zi,w_out)
      end