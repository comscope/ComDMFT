      subroutine prnatom(key,isort,nsp,acc,dff,nds)
c	key=0  - atom
c	key=1  - solid
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: key,isort,nsp,nds(maxntle,nrel*maxb+1,nsp)
      real*8, intent(in)  :: acc(maxntle,nrel*maxb+1,nsp),
     &                       dff(maxntle,nrel*maxb+1,nsp)
	integer :: isp,li,ival,l,ntl,ie,mt,mtc,ir,itc,jtc,i_len
      real*8 :: c_ov,c2,dqdall
	real*8, allocatable :: work(:)
	if(.not.maswrk) return
	allocate(work(0:maxnrad))
      c2=clight*clight
	do isp=1,nsp
	  write(iun,*)' Valence states'
	  if(key.eq.0) then
	    write(iun,201)
201	    format(' li',6x,'E',14x,'Psi_MT')
          ival=0
	    do li=1,nrel*lmb(isort)+1
	      if(irel.le.1) then
		    l=li-1
	      else
	        l=li/2
	      endif
            ntl=ntle(l,isort)
            do ie=1,ntl
	        ival=ival+1
	        write(iun,202)li,eny(ie,li,isort,isp)
	      enddo
202	      format(i3,2g14.7)
	    enddo
c ---------------------------------------------------------------------------------
	  else if(key.eq.1) then
c ---------------------------------------------------------------------------------
	    write(iun,*)' Potential parameters'
	    write(iun,103)
103	    format(4x,' li',4x,' ie',5x,'accuracy',8x,'Dnu',
     &           12x,'Enu',7x,' Norm')
	    do li=1,nrel*lmb(isort)+1
	      if(irel.le.1) then
		    l=li-1
	      else
	        l=li/2
	      endif
            do ie=1,ntle(l,isort)
	        write(iun,104)li,ie,acc(ie,li,isp)-1.d0,
     &			            dny(ie,li,isort,isp),eny(ie,li,isort,isp),
     &                      dff(ie,li,isp)
104 	        format(2(4x,i3),1x,e12.5,3(1x,g14.7))
	      enddo
	    enddo
	    write(iun,105)
105	    format(4x,' li',4x,' ie',5x,'fi',10x,'dfi',12x,'Core_ovrl',5x,
     &           'Nodes')
	    do li=1,nrel*lmb(isort)+1
	      if(irel.le.1) then
		    l=li-1
	      else
	        l=li/2
	      endif
            ntl=ntle(l,isort)
	      mtc=ntle_cor(li,isort)  !! last core level for this l
            do ie=1,ntl
	        mt=indfun0(1,ie,li,isort)
	        mt=ind_wf(mt,isort)
	        if(mtc==0) then
	          c_ov=0.d0
	        else
	          itc=indcor0(mtc,li,isort)
	          jtc=indcor(itc,isp,isort)
			    do ir=0,nrad(isort)
	            work(ir)=pcor(jtc+ir)*gfun(mt+ir,isp)
	            if(irel.ge.1) work(ir)=work(ir)
     &	              +qcor(jtc+ir)*gfund(mt+ir,isp)/c2
	                work(ir)=work(ir)*r(ir,isort)*dr(ir,isort)
	          enddo
                c_ov=dqdall(h(isort),work,nrad(isort))
	        endif
	        write(iun,106)li,ie,fi(ie,li,isort,isp),
     &		                dfi(ie,li,isort,isp),c_ov,nds(ie,li,isp)
106 	        format(2(4x,i3),3(1x,e12.5),7x,i3)
c ---------- Plotting the correlated orbital --------------------------
              if(correlated(ie,l,isort)=='L'.or.
     &           correlated(ie,l,isort)=='A') then
	          if(maswrk) then
	            i_len=len_trim(allfile)
	            if(isort==1) open(3,file=allfile(1:i_len)//'orb1.rad')
	            if(isort==2) open(3,file=allfile(1:i_len)//'orb2.rad')
	            if(isort==3) open(3,file=allfile(1:i_len)//'orb3.rad')
	            do ir=1,nrad(isort)
	              write(3,*)r(ir,isort),p_f(ir,ie,li,isort,isp)
     &	                                 /r(ir,isort)
	            enddo
		          close(3)
	          endif
	        endif
	      enddo  !! over ie
	    enddo
	  endif
      enddo   !! over isp
	deallocate(work)
      end