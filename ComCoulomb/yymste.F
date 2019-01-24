      subroutine yymste(b1,lmaxa,lmaxb,nlma0,nlmb0,b,nla,nlb,ila0,ilb0,
     &  eiphi)

c     multiplies by phase and copies into arrays at ila0,ilb0.

      implicit none
      integer :: lmaxa,lmaxb,nlma0,nlmb0,nla,nlb,ila0,ilb0,i,j,nlma,nlmb
      complex*16 :: b1(nlma0,nlmb0),eiphi,b(nla,nlb)

      nlma=(lmaxa+1)**2
      nlmb=(lmaxb+1)**2
      do i=1,nlmb
        do j=1,nlma
          b(j+ila0,i+ilb0)=eiphi*b1(j,i)
        enddo
      enddo

      end
