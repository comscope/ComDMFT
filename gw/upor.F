      subroutine upor(n,a)
	implicit none
	integer, intent(in) :: n
	real*8, intent(out) :: a(n)
	integer :: i,k,j
      real*8 :: s
      do i=1,n-1
         k=i
         s=a(i)
         do j=i+1,n
            if(a(j).lt.s) then
               s=a(j)
               k=j
            endif
         enddo
         if(k.ne.i) then
            a(k)=a(i)
            a(i)=s
         endif
      enddo
      end
      
      
      subroutine upor_c(n,n0,e,vl,vr)
	implicit none
	integer :: i,k,j,n,n0
	real*8 :: a(n),s
      complex*16 :: e(n),vl(n0,n0),vr(n0,n0),tmp(n),t
	do i=1,n
	  a(i)=dreal(e(i))
	enddo
      do i=1,n-1
         k=i
         s=a(i)
         do j=i+1,n
            if(a(j).lt.s) then
               s=a(j)
               k=j
            endif
         enddo
         if(k.ne.i) then
            a(k)=a(i)
            a(i)=s
	      t=e(k)
	      e(k)=e(i)
	      e(i)=t
	      call zcopy(n,vl(1,k),1,tmp,1)
	      call zcopy(n,vl(1,i),1,vl(1,k),1)
	      call zcopy(n,tmp,1,vl(1,i),1)
	      call zcopy(n,vr(1,k),1,tmp,1)
	      call zcopy(n,vr(1,i),1,vr(1,k),1)
	      call zcopy(n,tmp,1,vr(1,i),1)
         endif
      enddo
      end
      
      
      subroutine degen_list(eps,n,a,ndeg,nn_deg,n01)
	implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: eps,a(n)
	integer, intent(out) :: ndeg,nn_deg(n),n01(2,n)
	integer :: i
	ndeg=1
	nn_deg(ndeg)=1
	n01(:,1)=1
	i=2
1	if(abs(a(i)-a(i-1))<eps) then
        nn_deg(ndeg)=nn_deg(ndeg)+1
        n01(2,ndeg)=i
      else
        ndeg=ndeg+1
        nn_deg(ndeg)=1
        n01(:,ndeg)=i
      endif
      i=i+1
      if(i<=n) goto 1
      end