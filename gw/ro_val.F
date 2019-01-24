      subroutine ro_val
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: ispin,isort,isym,mt,km,km1,iwf1,ir,isymb,mtb,iwf,isp,
     &	         jsort,mt1
	real*8 :: c2,wk,tmp
	c2=clight**2
	do ispin=1,nspin
        isp=min(ispin,nspin_1)
	  do isort=1,nsort
	    do isym=1,nsym(isort)
	      mt=indmt(isym,isort,ispin)
            do km=1,lfun(isort)
	        iwf=ind_wf(km,isort)
              do km1=1,lfun(isort)
	          iwf1=ind_wf(km1,isort)
	          if(abs(t_ro(km,km1,isym,isort,ispin)).gt.1.d-10) then
                  do ir=0,nrad(isort)
                    wk=gfun(iwf+ir,isp)*gfun(iwf1+ir,isp)
                    ro(mt+ir)=ro(mt+ir)+t_ro(km,km1,isym,isort,ispin)*wk
	            enddo
	            if(irel.ge.1) then
                    do ir=0,nrad(isort)
                      wk=gfund(iwf+ir,isp)*gfund(iwf1+ir,isp)/c2
                      ro(mt+ir)=ro(mt+ir)
     &					       +t_ro(km,km1,isym,isort,ispin)*wk
	              enddo
	            endif
	          endif
	        enddo   !! over km1
	      enddo   !!  over km
	    enddo   !! over isym
	    if(magn.eq.2) then
	      do isymb=1,nsymb(isort)
	        mtb=indmtb(isymb,isort)
              do km=1,lfun(isort)
	          iwf=ind_wf(km,isort)
                do km1=1,lfun(isort)
	            iwf1=ind_wf(km1,isort)
	            if(abs(t_m(km,km1,1,isymb,isort)).gt.1.d-10) then
                    do ir=0,nrad(nsort)
                      wk=gfun(iwf+ir,isp)*gfun(iwf1+ir,isp)
                      spmt(mtb+ir)=spmt(mtb+ir)
     &				          +t_m(km,km1,1,isymb,isort)*wk
	              enddo
	            endif
	            if(abs(t_m(km,km1,2,isymb,isort)).gt.1.d-10) then
                    do ir=0,nrad(nsort)
                      wk=gfund(iwf+ir,isp)*gfund(iwf1+ir,isp)
                      spmt(mtb+ir)=spmt(mtb+ir)
     &				          +t_m(km,km1,2,isymb,isort)*wk
	              enddo
	            endif
	          enddo   !! over km1
	        enddo   !!  over km
	      enddo   !! over isymb
	    endif
	  enddo  !! over isort
	enddo   !! over isp
! --- IREL<2 AFM symmetrization --------------------------------
      if(irel<2.and.nspin==2) then
        do isort=1,nsort
          if(iafm(isort)==0) cycle
          jsort=iafm(isort)
	    mt=indmt(1,isort,1)
	    mt1=indmt(1,jsort,2)
          do ir=0,nrad(isort)
            tmp=0.5d0*(ro(mt+ir)+ro(mt1+ir))
            ro(mt+ir)=tmp
            ro(mt1+ir)=tmp
	    enddo
	  enddo
	endif
      end