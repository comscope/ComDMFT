      subroutine fgtlind
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
	integer :: isort,ind,ndimv,lm,i
      maxfgtl=0
      do isort=1,nsort
        ind=0
        ndimv=(lmpb(isort)+1)**2
        do lm=1,ndimv
          if(sym(lm,isort)) then
            ind=ind+1
            indfgtl(lm,isort)=ind
            goto 1
          endif
          if(magn.eq.2) then
            do i=1,3
              if(symb(lm,i,isort)) then
                ind=ind+1
                indfgtl(lm,isort)=ind
                goto 1
              endif
            enddo
          endif
          indfgtl(lm,isort)=0
1         continue
        enddo  !!! over lm
        maxfgtl=max(maxfgtl,ind)
      enddo   !!! over iatom
      end