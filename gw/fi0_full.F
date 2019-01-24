      subroutine fi0_full(fifim,isp1,isp2)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: isp1,isp2
      real*8, intent(out) :: fifim(maxel,maxel,n_pbmtm,nsort)
	integer :: isort,ndimb,ndimpb,lm3,l3,nm3,lm,l,nm,ii,iii,indx,lget,
     &	         icg1,icg2,icg,lmi,l_i,iei,lfi,ie3,km3,lf3,ie,km,lf,
     &           in3,in,jj3,jj,i,li3,ii3,mj3,li,mj
      real*8 :: cgg
	fifim=0.d0
	do isort=1,nsort
	  ndimb=nrel*(lmb(isort)+1)**2
	  ndimpb=(lmpb(isort)+1)**2
	  do lm3=1,ndimb
	    if(irel.le.1) then
	      l3=lget(lm3)
	      li3=l3+1
	    else if(irel.eq.2) then
	      call getlimj(lm3,l3,ii3,mj3,li3,0)
	    endif
	    nm3=ntle(l3,isort)
	    do lm=1,ndimb
	      if(irel.le.1) then
	        l=lget(lm)
	        li=l+1
	      else if(irel.eq.2) then
	        call getlimj(lm,l,ii,mj,li,0)
	      endif
	      nm=ntle(l,isort)
	      if(irel.le.1) then
		      ii = max0(lm3,lm)
		      iii = min0(lm3,lm)
              indx = (ii*(ii-1))/2 + iii
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1) - 1 
	      else
		      icg1=indxcgr(lm3,lm,1,0)
              icg2=indxcgr(lm3,lm,2,0)
	      endif
            do icg = icg1, icg2
              if(irel.le.1) then
		        lmi = jcg(icg)
	          cgg=cg(icg)
	        else
	          lmi=jcgr(icg)
	          cgg=cgr(icg)
	        endif
	        if(lmi.le.ndimpb) then
	          l_i=lget(lmi)
	          do iei=1,ntle_pb(l_i,isort)
	            lfi=ind_prod0(iei,l_i,isort)
	            i=ind_prod(iei,lmi,isort)
			      do ie3=1,nm3
			        in3=1
			        if(augm(ie3,l3,isort)/='LOC') in3=2
			        do jj3=1,in3
	                km3=indbasa(jj3,ie3,lm3,isort)
	                lf3=lf_isz(km3,isort)
			          do ie=1,nm
			            in=1
			            if(augm(ie,l,isort)/='LOC') in=2
			            do jj=1,in
	                    km=indbasa(jj,ie,lm,isort)
	                    lf=lf_isz(km,isort)
	                    fifim(km,km3,i,isort)=
     &                          cgg*fi0(lf3,lf,lfi,isort,isp1,isp2)
	                  enddo
	                enddo   !! over ie
	              enddo
	            enddo   !! over ie3
	          enddo   !!! over iei
	        endif
	      enddo   !! over icg
	    enddo   !! over lm
	  enddo  !! over lm3
	enddo  !! over isort
      end
