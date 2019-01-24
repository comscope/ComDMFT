      subroutine vcoulat(z,a,b,rofi,rho,nr,vhrmax,v)
c  hartree potential for spherical rho. vhrmax=value at rmax.
c  returns vsum = integral over pot.
c  rho = 'spherical chden' = 4pi*rhotrue. v=vtrue.
	use units_mod
      implicit none
	integer, intent(in) :: nr
      real*8, intent(in) :: rofi(0:nr),rho(0:nr),z,a,b,vhrmax
      real*8, intent(out) :: v(0:nr)
	integer :: ir
	real*8 :: sqpi4,rmax,r2,r3,r4,f2,f3,f4,x23,x34,cc,bb,dd,a2b4,df,r,
     &	      drdi,srdrdi,g,f,y2,y3,ro,y4,vnow
	sqpi4=sqrt(4.d0*pi)
	rmax=rofi(nr)
c ----- approximate rho by cc+bb*r+dd*r*r near zero ---
      r2=rofi(1)
      r3=rofi(2)
      r4=rofi(3)
      f2=rho(1)
      f3=rho(2)
      f4=rho(3)
      x23=(r3*r3*f2-r2*r2*f3)/(r3-r2)
      x34=(r4*r4*f3-r3*r3*f4)/(r4-r3)
      cc=(r2*x34-r4*x23) / (r3*(r2-r4))
      bb=((r2+r3)*x34-(r3+r4)*x23) / (r3*r3*(r4-r2))
      dd=(f2-bb*r2-cc)/r2**2
c ------ numerov for inhom solution --------
      a2b4=a*a/4.d0
      v(0)=1.d0
      df=0.d0
      do ir=1,2
        r=rofi(ir)
        drdi=a*(r+b)
        srdrdi=dsqrt(drdi)
        v(ir)=v(0)-r*r*(cc/3.d0+r*bb/6.d0+r*r*dd/10.d0)
        g=v(ir)*r/srdrdi
        f=g*(1.d0-a2b4/12.d0)
        if(ir.eq.1) y2=-2.d0*f2*r2*drdi*srdrdi
        if(ir.eq.2) y3=-2.d0*f3*r3*drdi*srdrdi
        df=f-df
      enddo
      ir=2
  13  ir=ir+1
      r=rofi(ir)
      drdi=a*(r+b)
      srdrdi=dsqrt(drdi)
      ro=rho(ir)
      y4=-2.d0*drdi*srdrdi*ro*r
      df=df+g*a2b4+(y4+10.d0*y3+y2)/12.d0
      f=f+df
      g=f/(1.d0-a2b4/12.d0)
      v(ir)=g*srdrdi/r
      y2=y3
      y3=y4
      if(ir.lt.nr) goto 13
	call dscal(nr+1,4.d0*pi,v(0),1)
c ------ add constant to get v(nr)=vhrmax -------
      vnow=v(nr)-2.d0*z/rmax*sqpi4
      do ir=0,nr
        v(ir)=v(ir)+(vhrmax-vnow)
      enddo
      end