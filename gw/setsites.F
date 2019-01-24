      subroutine setsites(delta,irs,ns,nvcmax,wd,nsitmax,ierr,radius)
***********************************************************************
*   Set nearest site distances for every kind of atom.                *
***********************************************************************
	use atom_mod
	use manager_mod
      use sitegen_mod
	use solid_mod
      implicit none
      integer, intent(in) :: nvcmax,nsitmax
      real*8, intent(in) :: radius
      integer, intent(out) :: irs(nsitmax,nsort),ns(nsort),ierr
	real*8, intent(out) :: wd(nsitmax,maxnsym,nsort),
     &	                   delta(nsitmax,nsort)
	integer :: isort,iatom,jatom,ic,isn,isite,key,jsn,lm,isym,nsite
	real*8 :: rmax,dels,delx,dely,delz
      integer, allocatable :: iref(:)
	real*8, allocatable :: ylm(:),site(:,:),taut(:,:)
        allocate(lastg_site(nvcmax/2))
       	allocate(iref_vgen(nvcmax))
	allocate(vect_site(3,0:nvcmax))
	allocate(rco_site(nvcmax/2))
	allocate(tauc_site(3,natom))
	allocate(vec_vgen(3,nvcmax))
	allocate(rw_vgen(nvcmax))
	allocate (ylm(limlpb),site(3,nvcmax),taut(3,natom),iref(nvcmax))
      ierr=0
	do isort=1,nsort
        do iatom=1,natom
          if(is(iatom).eq.isort) goto 1
	  enddo
1	  continue
c  calculate all vectors containing in a sphere surrounding given atom
        rmax=radius/par
        do jatom=1,natom
          do ic=1,3
            taut(ic,jatom)=tau(ic,jatom)-tau(ic,iatom)
          enddo
        enddo
        call sitegen(rmax,rbas(1,1),rbas(1,2),rbas(1,3),taut,natom,
     &	           nsite,site,iref,nvcmax)
c  sort vectors according to their length
        isn=0
        do isite=1,nsite
          dels=sqrt(site(1,isite)**2+site(2,isite)**2+
     &              site(3,isite)**2)
          if(dels.eq.0.d0) goto 50
          delx=site(1,isite)/dels
          dely=site(2,isite)/dels
          delz=site(3,isite)/dels
c  find spherical harmonics (wigner's matrices when m=0)
          call sphharm(delx,dely,delz,lmpb(isort),ylm)
          key=1
c  compare with previous vectors
          do jsn=1,isn
            if(delta(jsn,isort).eq.dels*par.and.
     &         irs(jsn,isort).eq.iref(isite)) then
              key=0
              do lm=1,(lmpb(isort)+1)**2
                if(sym(lm,isort)) then
                  isym=lmsym(lm,isort)
                  wd(jsn,isym,isort)=wd(jsn,isym,isort)+ylm(lm)
                endif
              enddo
              goto 50
            endif
          enddo
          if(key.ne.1) goto 50
c  if none add to array
          isn=isn+1
          if(isn.gt.nsitmax) ierr=1
          if(isn.gt.nsitmax) return
          delta(isn,isort)=dels*par
          do lm=1,(lmpb(isort)+1)**2
            if(sym(lm,isort)) then
              isym=lmsym(lm,isort)
              wd(isn,isym,isort)=ylm(lm)
            endif
          enddo
          irs(isn,isort)=iref(isite)
50        continue
        enddo     !!! over isite
        ns(isort)=isn
      enddo    !!! over isort
	deallocate (ylm,taut,site,iref)
        deallocate(lastg_site,vect_site,rco_site,tauc_site,iref_vgen,
     &             vec_vgen,rw_vgen)
      end
