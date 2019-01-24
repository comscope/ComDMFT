      subroutine zx_hf(fi0m,zz,zx,zxx,ispin)
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
	integer, intent(in) :: ispin
	integer :: iatom,isort,nl,ind_lm,ind_lme,ndimpb,i,ind_klm,ndimb,
     &           lm1,lm2,l1,lget,l2,ii,iii,indx,icg1,icg2,i2,mj2,li2,
     &           i1,li1,icg,lmi,mj1,l_i,nm1,nm2,iei,lfi,km,ie2,km2,
     &           lf2,ie1,km1,lf1,in,jn,in1,jn1
      complex*16 :: zz(nfun),fi0m(maxel,maxel,nlb2),zx(maxel,nlb2),
     &	            zxx(maxel,n_pbmt),cgg,cc
	zxx=(0.d0,0.d0)
	do iatom=1,natom
        isort=is(iatom)
        nl=lfunm(isort)
	  ind_lme=io_lem(iatom)
c ------------ ZX -----------------------------------------------
	  ind_lm=iob2(iatom)-1
	  ndimpb=(lmpb(isort)+1)**2
	  do i=ind_lm+1,ind_lm+ndimpb
	    call zgemv('n',nl,nl,(1.d0,0.d0),fi0m(1,1,i),maxel,
     &		         zz(ind_lme),1,(0.d0,0.d0),zx(1,i),1)
	  enddo   !! over i
c ------------ ZXX ----------------------------------------------
	  ind_klm=iopb(iatom)-1
	  ndimb=nrel*(lmb(isort)+1)**2
	  do lm1=1,ndimb
	    do lm2=1,ndimb
	      if(irel<=1) then
	        l1=lget(lm1)
	        l2=lget(lm2)
		      ii = max0(lm1,lm2)
		      iii = min0(lm1,lm2)
              indx = (ii*(ii-1))/2 + iii
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1) - 1 
	      else
	        call getlimj(lm2,l2,i2,mj2,li2,0)
	        call getlimj(lm1,l1,i1,mj1,li1,0)
		      icg1=indxcgr(lm1,lm2,1,0)
              icg2=indxcgr(lm1,lm2,2,0)
	      endif
            do icg = icg1, icg2
              if(irel<=1) then
		        lmi = jcg(icg)
	          cgg=dcmplx(cg(icg),0.d0)
	        else
	          lmi=jcgr(icg)
	          cgg=cgr(icg)
	        endif
	        if(lmi.le.ndimpb) then
	          l_i=lget(lmi)
	          nm1=ntle(l1,isort)
	          nm2=ntle(l2,isort)
	          do iei=1,ntle_pb(l_i,isort)
	            lfi=ind_prod0(iei,l_i,isort)
	            km=ind_klm+ind_prod(iei,lmi,isort)
			      do ie2=1,nm2
			        in=1
			        if(augm(ie2,l2,isort)/='LOC') in=2
			        do jn=1,in
	                km2=indbasa(jn,ie2,lm2,isort)
	                lf2=lf_isz(km2,isort)
	                cc=cgg*zz(ind_lme+km2-1)
	                do ie1=1,nm1
			            in1=1
			            if(augm(ie1,l1,isort)/='LOC') in1=2
			            do jn1=1,in1
	                    km1=indbasa(jn1,ie1,lm1,isort)
	                    lf1=lf_isz(km1,isort)
	                    zxx(km1,km)=zxx(km1,km)
     &	                             +cc
     &                               *fi0(lf1,lf2,lfi,isort,ispin,ispin)
                        enddo
	                enddo  ! over ie1
	              enddo
	            enddo  ! over ie2
	          enddo  !  over iei
	        endif
	      enddo   ! over icg
	    enddo   ! over lm2
	  enddo    ! over lm1
	enddo   !! over iatom
      end