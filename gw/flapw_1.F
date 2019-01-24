      subroutine flapw_1
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
      integer :: isort,i,lm,iatom,nbn,ndimb,li,l,ndim,m,ie,lf,jj,ii,
     &	         iend,mj,km,mt,icor,ilm,n,j,ispin,iwf,jn
      do isort=1,nsort
c --------- First atoms of the classes ----------------------------------
        do iatom=1,natom
	    if(is(iatom).eq.isort) then
	      iat_1(isort)=iatom
	      exit
	    endif
	  enddo
      enddo
      do isort=1,nsort
        i=0
        lm=0
        nbn=nrel*lmb(isort)+1
        ndimb=nrel*(lmb(isort)+1)**2
        do li=1,nbn
          if(irel.ne.2) l=li-1
          if(irel.eq.2) l=li/2
          if(irel.ne.2) ndim=li+li-1
          if(irel.eq.2) ndim=((li+1)/2)*2
          do m=1,ndim
            lm=lm+1
            do ie=1,ntle(l,isort)
              if(augm(ie,l,isort)=='LOC') then
                i=i+1
                indbasa(1,ie,lm,isort)=i
              endif
            enddo
          enddo
        enddo
        do jn=1,2
          lm=0
          do li=1,nbn
            if(irel.ne.2) l=li-1
            if(irel.eq.2) l=li/2
            if(irel.ne.2) ndim=li+li-1
            if(irel.eq.2) ndim=((li+1)/2)*2
            do m=1,ndim
              lm=lm+1
              do ie=1,ntle(l,isort)
                if(augm(ie,l,isort)=='MTO') then
                  i=i+1
                  indbasa(jn,ie,lm,isort)=i
                endif
              enddo
            enddo
          enddo
        enddo
        do jn=1,2
          lm=0
          do li=1,nbn
            if(irel.ne.2) l=li-1
            if(irel.eq.2) l=li/2
            if(irel.ne.2) ndim=li+li-1
            if(irel.eq.2) ndim=((li+1)/2)*2
            do m=1,ndim
              lm=lm+1
              do ie=1,ntle(l,isort)
                if(augm(ie,l,isort)=='APW') then
                  i=i+1
                  indbasa(jn,ie,lm,isort)=i
                endif
              enddo
            enddo
          enddo
        enddo
c ------------------------------------------------------------        
	  lf=0
        do li=1,nrel*lmb(isort)+1
	    if(irel.lt.2) then
	      l=li-1
	      jj=l
	    else
	      l=li/2
	      ii=mod(li,2)
	      if(ii.eq.0) ii=-1
            jj=l+l+ii
	    endif
          do ie=1,ntle(l,isort)
            iend=1
            if(augm(ie,l,isort)/='LOC') iend=2
            do i=1,iend
	        lf=lf+1
	        ie_val(lf,isort)=ie
	        li_val(lf,isort)=li
	        in_val(lf,isort)=i
              do mj=-jj,jj,nrel
                if(irel.lt.2) then
			      lm=l*(l+1)+mj+1
	          else
		          call getlimj(lm,l,ii,mj,li,1)
	          endif
                km=indbasa(i,ie,lm,isort)
	          lm_isz(km,isort)=lm
	          lf_isz(km,isort)=lf
	          ie_isz(km,isort)=ie
	        enddo
	      enddo
	    enddo
	  enddo
	enddo
c -------  Core information ----------------------------------
	ntle_cor=0
	mt=1
      do isort=1,nsort
	  icor=0
	  ilm=0
	  do li=1,nrel_core*lmb(isort)+1
          if(irel_core.ne.2) l=li-1
          if(irel_core.eq.2) l=li/2
	    do n=l+1,konfig(l,isort)-1
	      ntle_cor(li,isort)=ntle_cor(li,isort)+1
	      icor=icor+1
	      indcor0(n-l,li,isort)=icor
	      li_cor(icor,isort)=li
            nod_cor(icor,isort)=n-l-1
	      if(irel_core.le.1) then
	        do m=-l,l
	          ilm=ilm+1
	          lm=l*(l+1)+m+1
			    lm_cor(ilm,isort)=lm
			    lf_cor(ilm,isort)=icor
	          indcore(n-l,lm,isort)=ilm
			  enddo
		    else if(irel_core.eq.2) then
	        ii=li-2*l
	        if(ii.eq.0) ii=-1
	        j=2*l+ii
	        do mj=-j,j,2
	          ilm=ilm+1
		        call getlimj(lm,l,ii,mj,li,1)
			    lm_cor(ilm,isort)=lm
			    lf_cor(ilm,isort)=icor
	          indcore(n-l,lm,isort)=ilm
			  enddo
		    endif
	    enddo  !! over n
	  enddo  !! over li
	  ii=nrad(isort)+1
        do ispin=1,nspin_0
	    do icor=1,ncor(isort)
            indcor(icor,ispin,isort)=mt
            mt=mt+ii
	    enddo
	  enddo
      enddo
	iwf=1
	do isort=1,nsort
        do i=1,lfun(isort)
	    ind_wf(i,isort)=iwf
	    iwf=iwf+nrad(isort)+1
	  enddo
	enddo
	maxwf=iwf-1
      end