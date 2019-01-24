      subroutine tetrahedron(et,z,ar)
********************************************************
*  Integration over one tetrahedron.                   *
*  See Ph.Lambin and J.P.Vigneron, PRB 29, 3430 (1984) *
********************************************************
      implicit none
      real*8 :: et(4),a,b
	complex*16 :: z,ar(4)
	logical :: eqv
	eqv(a,b)=(abs(a-b).lt.1.d-12)
c -----------------------------------------------------------------------
	if(eqv(et(1),et(2))) then
	  if(eqv(et(2),et(3))) then
	    if(eqv(et(3),et(4))) then
	      call form_a6(z,et,ar)
	    else
	      call form_a3(z,et,ar)
	    endif
	  else
	    if(eqv(et(3),et(4))) then
	      call form_a4(z,et,ar)
	    else
	      call form_a12(3,4,1,2,z,et,ar)
	    endif
	  endif
	else
	  if(eqv(et(2),et(3))) then
	    if(eqv(et(3),et(4))) then
	      call form_a5(z,et,ar)
	    else
	      call form_a12(1,4,2,3,z,et,ar)
	    endif
	  else
	    if(eqv(et(3),et(4))) then
	      call form_a12(1,2,3,4,z,et,ar)
	    else
	      call form_7(z,et,ar)
	    endif
	  endif
	endif
      end
	subroutine form_a12(i,j,l,m,z,e,ar)
	implicit none
	integer :: i,j,l,m
	real*8 :: e(4)
	complex*16 :: z,ar(4),c1
c
	c1=2.d0*(z-e(m))/(e(i)-e(m))+(z-e(j))/(e(i)-e(j))
	c1=(1.d0,0.d0)+c1*log(z-e(i))
	ar(i)=(z-e(i))**2/(e(i)-e(m))**2/(e(j)-e(i))*c1
	c1=2.d0*(z-e(i))/(e(m)-e(i))+(z-e(j))/(e(m)-e(j))
	c1=(1.d0,0.d0)+c1*log(z-e(m))
	ar(i)=ar(i)+(z-e(m))**2/(e(m)-e(i))**2/(e(j)-e(m))*c1
	ar(i)=ar(i)+(z-e(j))**3/(e(j)-e(i))**2/(e(j)-e(m))**2*log(z-e(j))
c
	c1=2.d0*(z-e(m))/(e(j)-e(m))+(z-e(i))/(e(j)-e(i))
	c1=(1.d0,0.d0)+c1*log(z-e(j))
	ar(j)=(z-e(j))**2/(e(j)-e(m))**2/(e(i)-e(j))*c1
	c1=2.d0*(z-e(j))/(e(m)-e(j))+(z-e(i))/(e(m)-e(i))
	c1=(1.d0,0.d0)+c1*log(z-e(m))
	ar(j)=ar(j)+(z-e(m))**2/(e(m)-e(j))**2/(e(i)-e(m))*c1
	ar(j)=ar(j)+(z-e(i))**3/(e(i)-e(j))**2/(e(i)-e(m))**2*log(z-e(i))
c
	ar(l)=(z-e(j))**3/(e(j)-e(i))/(e(j)-e(m))**3*log(z-e(j))
	ar(l)=ar(l)+(z-e(i))**3/(e(i)-e(j))/(e(i)-e(m))**3*log(z-e(i))
	c1=((z-e(j))/(e(m)-e(j)))**2+((z-e(i))/(e(m)-e(i)))**2
	c1=c1+(z-e(j))*(z-e(i))/(e(m)-e(j))/(e(m)-e(i))
	c1=c1*log(z-e(m))+(0.5d0,0.d0)+(z-e(j))/(e(m)-e(j))
     &	                          +(z-e(i))/(e(m)-e(i))
	ar(l)=ar(l)+(z-e(m))/(e(m)-e(j))/(e(m)-e(i))*c1
	ar(m)=ar(l)
	end
	subroutine form_7(z,e,ar)
	implicit none
	integer :: i,k,j
	real*8 :: e(4),dme(4)
	complex*16 :: z,ar(4),ze(4)
	do i=1,4
	  ze(i)=z-e(i)
	  dme(i)=1.d0
	  do k=1,4
	    if(k.ne.i) dme(i)=dme(i)*(e(k)-e(i))
	  enddo
	enddo
	do i=1,4
	  ar(i)=ze(i)**2/dme(i)
	  do j=1,4
	    if(j.ne.i) then
	      ar(i)=ar(i)+ze(j)**3/dme(j)*log(ze(j)/ze(i))/(e(i)-e(j))
	    endif
	  enddo
	enddo
	end
	subroutine form_a3(z,e,ar)
	implicit none
	real*8 :: e(4),d43
	complex*16 :: z,ar(4),ze3,ze4
c
	d43=e(4)-e(3)
	ze3=z-e(3)
	ze4=z-e(4)
	ar(1)=ze4**3/d43**4*log(ze4/ze3)
	ar(1)=ar(1)+(6.d0*ze4**2-3.d0*d43*ze4+2.d0*d43**2)/6.d0/d43**3
	ar(2)=ar(1)
	ar(3)=ar(1)
	ar(4)=3.d0*ze4**2*ze3/d43**4*log(ze3/ze4)
	ar(4)=ar(4)-1.5d0*ze3*(2.d0*ze4-d43)/d43**3-1.d0/d43
	end
	subroutine form_a4(z,e,ar)
	implicit none
	real*8 :: e(4),d32
	complex*16 :: z,ar(4),ze3,ze2
c
	d32=e(3)-e(2)
	ze3=z-e(3)
	ze2=z-e(2)
	ar(1)=3.d0*ze3**2*ze2/d32**4*log(ze2/ze3)
	ar(1)=ar(1)-1.5d0*ze2*(2.d0*ze3-d32)/d32**3-1.d0/d32
	ar(2)=ar(1)
	ar(3)=3.d0*ze2**2*ze3/d32**4*log(ze3/ze2)
	ar(3)=ar(3)+1.5d0*ze3*(2.d0*ze2+d32)/d32**3+1.d0/d32
	ar(4)=ar(3)
	end
	subroutine form_a5(z,e,ar)
	implicit none
	real*8 :: e(4),d21
	complex*16 :: z,ar(4),ze1,ze2
c
	d21=e(2)-e(1)
	ze1=z-e(1)
	ze2=z-e(2)
	ar(1)=3.d0*ze1**2*ze2/d21**4*log(ze2/ze1)
	ar(1)=ar(1)+1.5d0*ze2*(2.d0*ze1+d21)/d21**3+1.d0/d21
	ar(2)=ze1**3/d21**4*log(ze1/ze2)
	ar(2)=ar(2)-(6.d0*ze1**2-3.d0*ze1*d21+2.d0*d21**2)/6.d0/d21**3
	ar(3)=ar(2)
	ar(4)=ar(2)
	end
	subroutine form_a6(z,e,ar)
	implicit none
	real*8 :: e(4)
	complex*16 :: z,ar(4)
c
	ar(1)=0.25d0/(z-e(1))
	ar(2)=ar(1)
	ar(3)=ar(1)
	ar(4)=ar(1)
	end