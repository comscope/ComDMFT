      subroutine memory_close
	use atom_mod
	use etot_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      if(allocated(gbs)) deallocate(gbs)
      if(allocated(igbs)) deallocate(igbs)
      if(allocated(rcosph)) deallocate(rcosph)
      if(allocated(lastg)) deallocate(lastg)
      if(allocated(tshift_r)) deallocate(tshift_r)
      end