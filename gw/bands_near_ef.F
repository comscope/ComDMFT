      subroutine bands_near_ef(iaa,ibb,nbn,ebn,mi,ma)
c     IAA > 0 - we take IBB valence!!! (i.e. after semicore) bands 
c               beginning from IAA
c     IAA = 0 - we take IBB bands closest to Ef
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: iaa,ibb,nbn
      real*8, intent(in) :: ebn(nbndf)
      integer, intent(out) :: mi,ma
	integer :: ib,imi,ima,first_bnd
	real*8 :: eps
	integer, allocatable :: iref(:)
	real*8, allocatable :: d(:)
	eps=1.d-6
	if(iaa>0) then
	  first_bnd=n_semic+iaa
        imi=first_bnd
        ima=first_bnd+ibb-1
        mi=imi
        do ib=imi-1,1,-1
          if(abs(ebn(imi)-ebn(ib))<eps) mi=ib
        enddo
        ma=ima
        do ib=ima+1,nbn
          if(abs(ebn(ima)-ebn(ib))<eps) ma=ib
        enddo	
      else if(iaa==0) then
	  allocate(d(nbn),iref(nbn))
	  do ib=1,nbn
	    d(ib)=abs(ebn(ib)-chem_pot)
	  enddo
		call ordering_all(nbn,d,iref)
c ---------- Check degeneracy --------------------------------
        imi=10000
        ima=-10000
        do ib=1,ibb
          imi=min(imi,iref(ib))
          ima=max(ima,iref(ib))
        enddo
        mi=imi
        do ib=imi-1,1,-1
          if(abs(ebn(imi)-ebn(ib))<eps) mi=ib
        enddo
        ma=ima
        do ib=ima+1,nbn
          if(abs(ebn(ima)-ebn(ib))<eps) ma=ib
        enddo
	  deallocate(d,iref)
      endif
      end