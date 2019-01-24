      subroutine zerov
***********************************************
*   set up muffin-tin zero for the potential  *
***********************************************
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: isort,mt,irad
	real*8 :: sqpi4,vmt0
      sqpi4=2.d0*sqrt(pi)
c  define mt-zero
      vmt0=0.d0
      do isort=1,nsort
        mt=indmt(1,isort,1)
        vmt0=vmt0+(v_mt_h(mt+nrad(isort))/sqpi4
     &           -2.d0*z(isort)/smt(isort))*nhsort(isort)
      enddo
      vmt0=vmt0/dfloat(natom)
c  shift potential according to new zero
      do isort=1,nsort
        mt=indmt(1,isort,1)
        do irad=mt,mt+nrad(isort)
          v_mt_h(irad)=v_mt_h(irad)-vmt0*sqpi4
        enddo
        vatmt(1,isort)=vatmt(1,isort)-vmt0*sqpi4
      enddo
      v_intr_h(1)=v_intr_h(1)-vmt0
      if(maswrk) write(iun,*)' MT-zero is ',vmt0
      end
