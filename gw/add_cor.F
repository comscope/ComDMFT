      subroutine add_cor(key)
c     key = -1 subtract core  
c     key =  1 add core      
	use atom_mod
	use etot_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: key
	integer :: isort,isp,mt
	do isort=1,nsort
	  do isp=1,nspin
          mt=indmt(1,isort,isp)
	    ro(mt:mt+nrad(isort))=ro(mt:mt+nrad(isort))
     &	                       +ro_core(0:nrad(isort),isort,isp)*key
	  enddo
	enddo  !! over isort
      end