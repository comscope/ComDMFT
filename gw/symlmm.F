      subroutine SYMLMM
****************************************************************
*  Define which coefficients in the expansion of the potential *
*  or charge density over spherical harmonics are not equal    *
*  to zero for a given crystalline group.                      *
*  First written by S.Yu. Savrasov                             *
****************************************************************
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer :: ind(17),isort,iatom,ndimv,lm,l,nls,m,m1,jatom,nls1,
     &  	       iwig,ig,mm,i,j
	real*8 :: ab
	real*8, allocatable :: tmp(:,:),tmp1(:,:,:)
      allocate(tmp(natom,-maxpb:maxpb),tmp1(natom,-maxpb:maxpb,3))
c   define which of the density or potential l-components will be zero
	do isort=1,nsort
        iatom=iat_1(isort)
        ndimv=(lmpb(isort)+1)**2
        nsym(isort)=0
        do lm=1,ndimv
          sym(lm,isort)=.false.
        enddo
        do l=0,lmpb(isort)
          nls=l*(2*l-1)*(2*l+1)/3+l+1
          do m=-l,l
            lm=l*(l+1)+m+1
            do m1=-l,l
              do jatom=1,natom
                tmp(jatom,m1)=0.d0
              enddo
            enddo
            nls1=nls+m
            do m1=-l,l
              iwig=nls1+(2*l+1)*(l+m1)
              do ig=1,ngroup
                jatom=ip(iatom,ig)
                tmp(jatom,m1)=tmp(jatom,m1)+u(iwig,ig)
              enddo
            enddo
            ab=0.d0
            do m1=-l,l
              do jatom=1,natom
                if(is(jatom).eq.isort) ab=ab+abs(tmp(jatom,m1))
              enddo
            enddo
            if(ab.gt.1.d-6) then
              nsym(isort)=nsym(isort)+1
              sym(lm,isort)=.true.
              lmsym(lm,isort)=nsym(isort)
            endif
          enddo   !!! over m
        enddo   !!! over l
      enddo   !!! over isort
c  print allowed lm-values
      maxnsym=0
      do isort=1,nsort
        if(maswrk.and.iprn.gt.0) then
          write(iun,*)txtel(isort)
          write(iun,1030)lmpb(isort)
          write(iun,1040)
        endif
        do l=0,min(lmpb(isort),8)
          mm=0
          do m=-l,l
            lm=l*(l+1)+m+1
            if(sym(lm,isort)) then
              mm=mm+1
              ind(mm)=m
            endif
          enddo
          if(mm.gt.0.and.maswrk.and.iprn.gt.0)
     &         write(iun,1050)l,(ind(i),i=1,mm)
        enddo
        if(maswrk.and.iprn.gt.0) write(iun,1070)nsym(isort)
        maxnsym=max(maxnsym,nsym(isort))
      enddo
      if(magn.eq.2) then
	  do isort=1,nsort
          iatom=iat_1(isort)
          ndimv=(lmpb(isort)+1)**2
          nsymb(isort)=0
          do i=1,3
            do lm=1,ndimv
              symb(lm,i,isort)=.false.
            enddo
          enddo
          do i=1,3
            do l=0,lmpb(isort)
              nls=l*(2*l-1)*(2*l+1)/3+l+1
              do m=-l,l
                lm=l*(l+1)+m+1
                do j=1,3
                  do m1=-l,l
                    do jatom=1,natom
                      tmp1(jatom,m1,j)=0.d0
                    enddo
                  enddo
                enddo
                nls1=nls+m
                do j=1,3
                  do m1=-l,l
                    iwig=nls1+(2*l+1)*(l+m1)
                    do ig=1,ngroup
                      jatom=ip(iatom,ig)
                      tmp1(jatom,m1,j)=tmp1(jatom,m1,j)+
     &                      u(iwig,ig)*ifia(ig)*agm(i,j,ig)
                    enddo
                  enddo
                enddo
                ab=0.d0
                do j=1,3
                  do m1=-l,l
                    do jatom=1,natom
                      if(is(jatom).eq.isort) ab=ab+abs(tmp1(jatom,m1,j))
                    enddo
                  enddo
                enddo
                if(ab.gt.1.d-6) then
                  nsymb(isort)=nsymb(isort)+1
                  symb(lm,i,isort)=.true.
                  lmsymb(lm,i,isort)=nsymb(isort)
                endif
              enddo   !!! over m
            enddo   !!! over l
          enddo   !!! over i
        enddo   !!! over isort
c  print allowed lm-values
        maxnsymb=0
        do isort=1,nsort
          if(maswrk.and.iprn.gt.0) write(iun,*)txtel(isort)
          do i=1,3
            if(i.eq.1.and.maswrk.and.iprn.gt.0) write(iun,1240)
            if(i.eq.2.and.maswrk.and.iprn.gt.0) write(iun,1250)
            if(i.eq.3.and.maswrk.and.iprn.gt.0) write(iun,1260)
            do l=0,min(lmpb(isort),8)
              mm=0
              do m=-l,l
                lm=l*(l+1)+m+1
                if(symb(lm,i,isort)) then
                  mm=mm+1
                  ind(mm)=m
                endif
              enddo
              if(mm.gt.0.and.maswrk.and.iprn.gt.0)
     &            write(iun,1050)l,(ind(j),j=1,mm)
            enddo
          enddo
          if(maswrk.and.iprn.gt.0) write(iun,1270)nsymb(isort)
          maxnsymb=max(maxnsymb,nsymb(isort))
        enddo
      endif   !!! only for magn=2
1030  format(1x,'Charge density is expanded in spherical harmonics ',
     &  'up to lmax=',i2)
1040  format(1x,'Non-zero elements of density allowed by symmetry are ',
     &            'the following:')
1050  format(1x,'l=',i2,' ; m=',17(1x,i3))
1070  format(1x,'Total # of non-zero components of density found ',i3)
1240  format(1x,'Non-zero elements of magnetic moment along X are ',
     &            'the following:')
1250  format(1x,'Non-zero elements of magnetic moment along Y are ',
     &            'the following:')
1260  format(1x,'Non-zero elements of magnetic moment along Z are ',
     &            'the following:')
1270  format(1x,'Total # of non-zero components of magn.mom. found ',i3)
	deallocate(tmp,tmp1)
      end