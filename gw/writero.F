      subroutine writero
************************************************************
*     store charge density at given iteration              *
************************************************************
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer :: i_len
      if(.not.maswrk) return
	i_len=len_trim(allfile)
      open(3,file=allfile(1:i_len)//'_el_density',form='unformatted')
c general information
      write(3)ro
	write(3)rointr
      if(magn.eq.2) then
        write(3)spmt
	  write(3)spintr
      endif
	close(3)
      end
      
      subroutine write_ro
************************************************************
*     store charge density at given iteration              *
************************************************************
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer :: i_len
      if(.not.maswrk) return
	i_len=len_trim(allfile)
      open(3,file=allfile(1:i_len)//'_ro',form='unformatted')
c general information
      write(3)ro
	write(3)rointr
      if(magn.eq.2) then
        write(3)spmt
	  write(3)spintr
      endif
	close(3)
      end
