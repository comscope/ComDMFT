      subroutine setk2(pnt_,gb,yr,yi,nbas,ind_y)
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
	integer, intent(in) :: nbas
      real*8, intent(in) :: pnt_(3),gb(3,nbas/nrel)
      integer, intent(out) :: ind_y(2*limlb,natom)
      real*8, intent(out) :: yr(*),yi(*)
	integer :: ind,iatom,isort,ndimb,lm
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      ind=0
      do iatom=1,natom
        isort=is(iatom)
        ndimb=nrel*(lmb(isort)+1)**2
        do lm=1,2*ndimb
          ind=ind+1
          ind_y(lm,iatom)=ind
        enddo
      enddo
      call gety(gb,yr,yi,pnt_,nbas,ind_y)
      end