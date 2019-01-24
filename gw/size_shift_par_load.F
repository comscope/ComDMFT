      subroutine size_shift_par_load(n,npr,ndim_pr,index_pr,ndim0)
      implicit none
	integer, intent(in) :: n,npr,ndim0
	integer, intent(out) :: ndim_pr(npr),index_pr(ndim0,npr)
	integer :: ipr,i
	ndim_pr=0
	do i=1,n
	  ipr=mod(i,npr)
	  if(ipr==0) ipr=npr
	  ndim_pr(ipr)=ndim_pr(ipr)+1
	  index_pr(ndim_pr(ipr),ipr)=i
	enddo
	end