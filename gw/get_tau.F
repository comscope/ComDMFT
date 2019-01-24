      subroutine get_tau
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
	integer :: iatom,i,j
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      do iatom=1,natom
        do i=1,3
          tau(i,iatom)=0.d0
          do j=1,3
            tau(i,iatom)=tau(i,iatom)+rbas(i,j)*stau(j,iatom)
          enddo
        enddo
      enddo
      end