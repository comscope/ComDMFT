      subroutine vexchat(rofi,rho,nr,v,rhomu,rhoeps,dr,h)
	use manager_mod
	use units_mod
      implicit none
	integer, intent(in) :: nr
      real*8, intent(in) :: rofi(0:nr),rho(0:nr),dr(0:nr)
      real*8, intent(inout) :: v(0:nr)
      real*8, intent(out) :: rhomu,rhoeps
	integer :: ir
	real*8 :: sqpi4,rho2,rho3,rho0,roup,rodn,gro,gmod,aamod,rolap,
     &	        EXCUP,EXCDN,VXCUP,VXCDN,dqdall,h
	real*8, allocatable :: work(:,:)
	allocate(work(0:nr,2))
	sqpi4=sqrt(4.d0*pi)
      rho2=rho(1)
      rho3=rho(2)
      rho0=(rho2*rofi(2)-rho3*rofi(1))/(rofi(2)-rofi(1))/sqpi4
      roup=rho0/2.d0
      rodn=rho0/2.d0
	call EXCVXC(ROUP,RODN,ilda,0,gro,gmod,aamod,rolap,
     &	          EXCUP,EXCDN,VXCUP,VXCDN)
	vxcup=vxcup*sqpi4
	excup=excup*sqpi4
      v(0)=v(0)+vxcup
	work(0,1)=vxcup*rho(0)*rofi(0)**2*dr(0)
	work(0,2)=excup*rho(0)*rofi(0)**2*dr(0)
      do ir=1,nr
        roup=rho(ir)/sqpi4/2.d0
        rodn=rho(ir)/sqpi4/2.d0
	  call EXCVXC(ROUP,RODN,ilda,0,gro,gmod,aamod,rolap,
     &	            EXCUP,EXCDN,VXCUP,VXCDN)
	  vxcup=vxcup*sqpi4
	  excup=excup*sqpi4
        v(ir)=v(ir)+vxcup
	  work(ir,1)=vxcup*rho(ir)*rofi(ir)**2*dr(ir)
	  work(ir,2)=excup*rho(ir)*rofi(ir)**2*dr(ir)
      enddo
      rhomu=dqdall(h,work(0,1),nr)
      rhoeps=dqdall(h,work(0,2),nr)
	deallocate(work)
      end