      subroutine xorder(nin,a,nout,b,iref)
**********************************************************
*       subroutine for ordering area with large          *
*           number of elements by decreasing             *
*                                                        *
*       input: nin  -number of elements,a(i=1,nin)-area. *
*      output: nout -number of non eqvivalent elements,  *
*              b(i=1,nout)-ordering area.                *
**********************************************************
      implicit none
      integer, intent(in) :: nin
      integer, intent(out) :: nout,iref(nin)
      real*8, intent(in) :: a(nin)
      real*8, intent(out) :: b(nin)
	integer :: j,nsum
	real*8 :: eps,s
      eps=1.d-8
c start up:
      iref=0
      nout=0
1     nout=nout+1
      s=1.d5
      do j=1,nin
         if(iref(j).eq.0) s=min(s,a(j))
      enddo
      b(nout)=s
      nsum=0
      do j=1,nin
         if(iref(j).eq.0.and.a(j)-s.lt.eps) iref(j)=nout
         if(iref(j).eq.0) nsum=nsum+1
      enddo
      if(nsum.gt.0) goto 1
      end