      subroutine inp_lda(irad,yl,rotot,isort,totm,sm)
	use atom_mod
	use etot_mod
	use manager_mod
	use solid_mod
      implicit none
      integer :: ndimv,ispin,lm,isym,mt,i,mtb,isort,irad,isymb
      real*8 :: yl(limlpb),rotot(2),totm(3),sm
	ndimv=(lmpb(isort)+1)**2
	if(nspin.eq.2) then
        do ispin=1,nspin
          rotot(ispin)=0.d0
          do lm=1,ndimv
            if(sym(lm,isort)) then
              isym=lmsym(lm,isort)
	        mt=indmt(isym,isort,ispin)
              rotot(ispin)=rotot(ispin)+ro(mt+irad)*yl(lm)
            endif
          enddo  !!! over lm
        enddo
	else if(nspin.eq.1) then
        rotot(1)=0.d0
        do lm=1,ndimv
          if(sym(lm,isort)) then
            isym=lmsym(lm,isort)
	      mt=indmt(isym,isort,1)
            rotot(1)=rotot(1)+ro(mt+irad)*yl(lm)/2.d0
          endif
        enddo  !!! over lm
	  if(magn.eq.1) then
c ------------ Copying to ispin=2 ------------------------------
          rotot(2)=rotot(1)
c -----------------------------------------------------------
        else if(magn.eq.2) then
c ----------- Calculate Full Magnetization -------------------
          do i=1,3
            totm(i)=0.d0
            do lm=1,ndimv
              if(symb(lm,i,isort)) then
                isymb=lmsymb(lm,i,isort)
	          mtb=indmtb(isymb,isort)
                totm(i)=totm(i)+spmt(mtb+irad)*yl(lm)
              endif
            enddo  !!! over lm
	    enddo
	    sm=sqrt(dot_product(totm,totm))
c ------- Forming RO_up and RO_down values ----------------------------
          rotot(2)=rotot(1)+0.5d0*sm
          rotot(1)=rotot(1)-0.5d0*sm
	  endif
	endif
      end