      SUBROUTINE EXCVXC2d(ROUP,RODN,ild0,EXCUP,EXCDN,VXCUP,VXCDN)
*************************************************************
*  Calculates 2D exchange-correlation energy and potential  *
*  after LSDA.                                              *
*    ROUP=Rho_up                                            *
*    RODN=Rho_dn                                            *
*  The ILDA=IXC-MOD(IXC,10) is avaiable with the following  *
*  parametrizations:                                        *
*  ILD0=0 - no exchange-correlation                         *
*      =1 - N.D.Drummond and R.J.Needs,                     *
*           Phys.Rev.B 79, 085414 (2009)
*      =6 - Original Gaspar-Kohn-Sham (rho**1/3)            *
*     In case ILD0<0 only correlation is taken              *
*************************************************************
      use etot_mod
      use manager_mod
      use units_mod
      IMPLICIT none
	integer, intent(in) :: ild0
      real*8, intent(in) :: roup,rodn
      real*8, intent(out) :: excup,excdn,vxcup,vxcdn
      integer :: ild
	real*8 :: exldaup,vxldaup,exldadn,vxldadn,rho,rs,eclda,vcldaup,
     & 	        vcldadn,a0,b0,c0,d0,e0,f0,g0,h0,aa,aad,ff,ffd,decdrs
C  Energies and potentials
      ild=iabs(ild0)
      EXCUP=0.D0
      EXCDN=0.D0
      VXCUP=0.D0
      VXCDN=0.D0
      IF(ild.EQ.0)RETURN
C  ********************************
C  LOCAL SPIN DENSITY APPROXIMATION
C  ********************************
      rho=roup+rodn
      rs=1.d0/sqrt(pi*rho)
C  Local exchange
      EXLDAUP=-4.d0*sqrt(2.d0)/3.d0/pi/rs
      VXLDAUP=-2.d0*sqrt(2.d0)/pi/rs
      EXLDAdn=-4.d0*sqrt(2.d0)/3.d0/pi/rs
      VXLDAdn=-2.d0*sqrt(2.d0)/pi/rs
******
      if(ild0>=0) then
        EXCUP=EXCUP+EXLDAUP
        EXCDN=EXCDN+EXLDADN
        VXCUP=VXCUP+VXLDAUP
        VXCDN=VXCDN+VXLDADN
      endif
C  Local correlation:
      ECLDA=0.D0
      VCLDAUP=0.D0
      VCLDADN=0.D0
******
      IF(ild==1) then
        a0=-0.1925d0
        b0=sqrt(2.d0)*(10.d0-3.d0*pi)/3.d0/pi
        c0=0.06979568d0
        d0=0.d0
        e0=1.0531003d0
        f0=0.04069122d0
        g0=0.3605953d0
        h0=0.d0
        aa=b0*rs+c0*rs*rs+d0*rs**3
        ff=e0*rs+f0*rs**1.5d0+g0*rs*rs+h0*rs**3
        eclda=a0+aa*log(1.d0+1.d0/ff)
        aad=b0+2.d0*c0*rs+3.d0*d0*rs**2
        ffd=e0+1.5d0*f0*rs**0.5d0+2.d0*g0*rs+3.d0*h0*rs**2
        decdrs=aad*log(1.d0+1.d0/ff)-aa/ff/(1.d0+ff)*ffd
        vcldaup=eclda-decdrs*rs/2.d0
        vcldadn=vcldaup
******
        EXCUP=EXCUP+ECLDA
        EXCDN=EXCDN+ECLDA
        VXCUP=VXCUP+VCLDAUP
        VXCDN=VXCDN+VCLDADN
      endif
      EXCUP=EXCUP*2.d0
      EXCDN=EXCDN*2.d0
      VXCUP=VXCUP*2.d0
      VXCDN=VXCDN*2.d0
      END