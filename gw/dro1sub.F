      subroutine dro1sub(drodr,d2rodr2,dspdr,d2spdr2,isort)
	use atom_mod
	use manager_mod
	use solid_mod
****************************************************************
*   calculate radial derivative of spherical part of rho       *
****************************************************************
      implicit none
	integer, intent(in) :: isort
      real*8, intent(out) :: drodr(0:nrad(isort),nsym(isort),nspin),
     &	                   d2rodr2(0:nrad(isort),nsym(isort),nspin),
     &                       dspdr(0:nrad(isort),*),
     &                       d2spdr2(0:nrad(isort),*)
	integer :: ndimv,ispin,lm,isym,mt,irad,i,isymb,mtb
	real*8 :: betta,hdif,diraddr
* define mesh parameters
      betta=r(nrad(isort),isort)/(exp(h(isort)*nrad(isort))-1.d0)
      hdif=1.d0/h(isort)
	ndimv=(lmpb(isort)+1)**2
      do ispin=1,nspin
        do lm=1,ndimv
          if(sym(lm,isort))then
            isym=lmsym(lm,isort)
* get first-order derivative
	      mt=indmt(isym,isort,ispin)
            call dif(ro(mt+1),nrad(isort)-1,drodr(1,isym,ispin))
            do irad=1,nrad(isort)
              diraddr=hdif/(r(irad,isort)+betta)
              drodr(irad,isym,ispin)=drodr(irad,isym,ispin)*diraddr
            enddo
* get second-order derivative
            call dif(drodr(1,isym,ispin),nrad(isort)-1,
     &               d2rodr2(1,isym,ispin))
            do irad=1,nrad(isort)
              diraddr=hdif/(r(irad,isort)+betta)
              d2rodr2(irad,isym,ispin)=d2rodr2(irad,isym,ispin)*
     &                                                       diraddr
            enddo
          endif
        enddo
      enddo
      if(magn.eq.2) then
	  do i=1,3
          do lm=1,ndimv
            if(symb(lm,i,isort))then
              isymb=lmsymb(lm,i,isort)
* get first-order derivative
	        mtb=indmtb(isymb,isort)
              call dif(spmt(mtb+1),nrad(isort)-1,dspdr(1,isymb))
              do irad=1,nrad(isort)
                diraddr=hdif/(r(irad,isort)+betta)
                dspdr(irad,isymb)=dspdr(irad,isymb)*diraddr
              enddo
* get second-order derivative
              call dif(dspdr(1,isymb),nrad(isort)-1,d2spdr2(1,isymb))
              do irad=1,nrad(isort)
                diraddr=hdif/(r(irad,isort)+betta)
                d2spdr2(irad,isymb)=d2spdr2(irad,isymb)*diraddr
              enddo
            endif
          enddo
	  enddo
      endif
      end