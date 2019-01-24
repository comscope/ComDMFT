      subroutine xxmste0(b1,b2,nlma0)

c     puts strux belonging to negative vector into b2,bb2

      implicit none
      integer :: nlma0,ilmb,lb,lget,ilma,la,sig
      complex*16 :: b1(nlma0,*),b2(nlma0,*)

      do ilmb=1,nlma0
        lb=lget(ilmb)
        do ilma=1,nlma0
          la=lget(ilma)
          sig=1-2*mod(la+lb,2)
          b2(ilma,ilmb)=sig*dconjg(b1(ilma,ilmb))
        enddo
      enddo
      end
