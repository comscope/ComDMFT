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
      open(3,file=
     $    allfile(1:i_len)//'_el_density',form='unformatted')
c general information
      write(3)ro
	write(3)rointr
      if(magn.eq.2) then
        write(3)spmt
	  write(3)spintr
      endif
	close(3)
      end
      
c$$$      subroutine write_ro
c$$$************************************************************
c$$$*     store charge density at given iteration              *
c$$$************************************************************
c$$$	use atom_mod
c$$$	use manager_mod
c$$$	use parallel_mod
c$$$	use solid_mod
c$$$      implicit none
c$$$      integer :: i_len
c$$$      if(.not.maswrk) return
c$$$	i_len=len_trim(allfile)
c$$$      open(3,file=allfile(1:i_len)//'_ro',form='unformatted')
c$$$c general information
c$$$      write(3)ro
c$$$	write(3)rointr
c$$$      if(magn.eq.2) then
c$$$        write(3)spmt
c$$$	  write(3)spintr
c$$$      endif
c$$$	close(3)
c$$$      end
