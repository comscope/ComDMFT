      subroutine get_nre(e,nre,l,rofi,v,z,nr)
      implicit none
	integer :: nre,nr,ir,l
      real*8 :: v(0:nr),rofi(0:nr),fllp1,alf,bet,e,z,ex
	fllp1=l*(l+1)
	do ir=nr,3,-1
	  alf=fllp1/rofi(ir)**2+v(ir)-2.d0*z/rofi(ir)-e
	  if(alf.le.0.d0) then
	    nre=ir
	    goto 1
	  endif
	  bet=sqrt(alf)
c	  if(bet*rofi(ir)<13.d0*log(10.d0)) then
c	    nre=ir
c	    goto 1
c	  endif
	  ex=exp(-bet*rofi(ir))
	  if(ex.gt.1.d-10) then
	    nre=ir
	    goto 1
	  endif
	enddo
1     continue
      nre=(nre/2)*2
      end