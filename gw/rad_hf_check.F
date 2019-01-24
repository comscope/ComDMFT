  	subroutine rad_hf_check(pn,qn,nre,node,nr)
      use manager_mod
      use units_mod
      implicit none
      integer, intent(in) :: nr
      integer, intent(inout) :: nre,node
      real*8, intent(inout) :: pn(0:nr),qn(0:nr)
	integer :: i,ii
      real*8 :: a
      integer, allocatable :: nn(:)
      real*8, allocatable :: p(:)
      allocate(p(node+1))
      allocate(nn(node+1))
      ii=1
      nn=0
      p=0.d0
      a=pn(1)
      do i=2,nre
        if(a*pn(i)<0.d0) then
          nn(ii)=i-1
          ii=ii+1
        endif
        a=pn(i)
        p(ii)=max(p(ii),abs(a))
      enddo
      nn(ii)=nre
      do i=2,ii
        if(p(i)/p(i-1)<1.d-1) then
          node=i-2
          nre=nn(i-1)
          goto 1
        endif
      enddo
1     deallocate(p,nn)
      pn(nre+1:nr)=0.d0
      qn(nre+1:nr)=0.d0
      if(abs(pn(nre))>1.d-7.and.nre<nr) nre=nre+1
	end