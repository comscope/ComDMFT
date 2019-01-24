      subroutine getind
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
	integer :: mt,ispin,isort,ii,isym
      allocate(indmt(maxnsym,nsort,nspin_0))
      mt=1
      maxmt=0
      do ispin=1,nspin_0
        do isort=1,nsort
          ii=nrad(isort)+1
          do isym=1,nsym(isort)
            indmt(isym,isort,ispin)=mt
            mt=mt+ii
            maxmt=maxmt+ii
          enddo
        enddo
      enddo
      if(magn.eq.2) then
        allocate(indmtb(maxnsymb,nsort))
        mt=1
        maxmtb=0
        do isort=1,nsort
          ii=nrad(isort)+1
          do isym=1,nsymb(isort)
            indmtb(isym,isort)=mt
            mt=mt+ii
            maxmtb=maxmtb+ii
          enddo
        enddo
      endif
      if(maxmtb.eq.0) maxmtb=1
      end