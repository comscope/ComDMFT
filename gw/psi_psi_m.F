      subroutine psi_psi_m(bbm,iq,n_pbt,ispin,n1,n)
c     Here K = Q
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: iq,n_pbt,ispin,n,n1
	complex*16, intent(out) :: bbm(n1,n,n_pbt)
      integer :: i,j,iatom,isort,ndimb,ndimpb,ind0,indp,lm1,l1,lget,i1,
     &           mj1,li1,nm1,lm,l,mj,li,nm,ii,iii,indx,icg1,icg2,icg,
     &           lmi,l_i,iei,lfi,ipb,ie1,km1,lf1,ie,km,lf,id,id0,ia0,
     &           ib0,ic0,jbas,jnd,ia,ib,ic,ibas,ind,ia1,ib1,ic1,in,jn,
     &           in1,jn1,irl,ibas0,iq0,isp
	complex*16 :: cgg,cc
	complex*16, allocatable :: znew(:,:),anew(:,:),s(:,:),tmp(:,:)
      isp=min(ispin,nspin_1)
      iq0=i_kref(iq)
	allocate(znew(nfun,n1))
	allocate(anew(nbasmpw,n1))
      call sym_z_0(znew,iq,z_bnd(1,1,iq0,ispin),n1,k_group(iq),
     &             pnt(1,iq))
	call sym_a_1(anew,iq,ev_bnd(1,1,iq0,ispin),n1,k_group(iq),iq0)
	allocate(s(nbndtm,nbndtm))
	allocate(tmp(nbndtm,nbndf))
	bbm=(0.d0,0.d0)
c ----------------- MT contribution ----------------------------------
	do iatom=1,natom
	  isort=is(iatom)
	  ndimb=nrel*(lmb(isort)+1)**2
	  ndimpb=(lmpb(isort)+1)**2
	  ind0=io_lem(iatom)-1   !! zero index in Z_BND
	  indp=iopb(iatom)-1   !! zero index in PB
	  do lm1=1,ndimb
	    if(irel.le.1) then
	      l1=lget(lm1)
	    else
	      call getlimj(lm1,l1,i1,mj1,li1,0)
	    endif
	    nm1=ntle(l1,isort)
	    do lm=1,ndimb
	      if(irel.le.1) then
	        l=lget(lm)
	      else
	        call getlimj(lm,l,i,mj,li,0)
	      endif
	      nm=ntle(l,isort)
	      if(irel.le.1) then
	        ii = max0(lm1,lm)
	        iii = min0(lm1,lm)
              indx = (ii*(ii-1))/2 + iii
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1) - 1 
	      else
		      icg1=indxcgr(lm,lm1,1,0)
              icg2=indxcgr(lm,lm1,2,0)
	      endif
            do icg = icg1, icg2
              if(irel.le.1) then
		        lmi = jcg(icg)
	          cgg=dcmplx(cg(icg),0.d0)
	        else
	          lmi=jcgr(icg)
	          cgg=cgr(icg)
	        endif
	        if(lmi.le.ndimpb) then
	          l_i=lget(lmi)
	          do iei=1,ntle_pb(l_i,isort)
	            lfi=ind_prod0(iei,l_i,isort)
	            ipb=indp+ind_prod(iei,lmi,isort)
		          do ie1=1,nm1
			        in1=1
			        if(augm(ie1,l1,isort)/='LOC') in1=2
			        do jn1=1,in1
	                km1=indbasa(jn1,ie1,lm1,isort)
	                lf1=lf_isz(km1,isort)
	                km1=ind0+km1
		              do ie=1,nm
			            in=1
			            if(augm(ie,l,isort)/='LOC') in=2
			            do jn=1,in
	                    km=indbasa(jn,ie,lm,isort)
	                    lf=lf_isz(km,isort)
	                    km=ind0+km
	                    cc=cgg*fi0(lf1,lf,lfi,isort,isp,isp)
	                    do j=1,n
	                      do i=1,n1
	                        bbm(i,j,ipb)=bbm(i,j,ipb)+
     &                                     cc*conjg(znew(km,i))
     &                                    *z_bnd(km1,j,1,ispin)
                            enddo
                          enddo
	                  enddo
	                enddo
	              enddo   !! over ie
	            enddo   !! over ie1
	          enddo   !!! over iei
	        endif
	      enddo   !! over icg
	    enddo   !! over lm
	  enddo  !! over lm1
	enddo  !! over iatom
c ----------------- Interstitial contribution ------------------------
	do id=1,nplwgw(iq0)
	  id0=indpw_gw(id,iq)
	  ipb=id+n_pbmt
	  ia0=igbs(1,id0)
	  ib0=igbs(2,id0)
	  ic0=igbs(3,id0)
c ------ We form the overlap matrix -----------------------------
	  do jbas=1,nbask(1)/nrel
          jnd=indgb(jbas,1)
          ia=ia0+igbs(1,jnd)
          ib=ib0+igbs(2,jnd)
          ic=ic0+igbs(3,jnd)
          do ibas=1,nbask(iq0)/nrel
            ind=indgb(ibas,iq)
            ia1=ia-igbs(1,ind)
            ib1=ib-igbs(2,ind)
            ic1=ic-igbs(3,ind)
            j=indplw(ia1,ib1,ic1)
	      if(complex_ro) then
	        s(ibas,jbas)=dcmplx(sovr(j),sovi(j))
	      else
	        s(ibas,jbas)=dcmplx(sovr(j),0.d0)
	      endif
	    enddo
	  enddo
	  do irl=1,nrel
          ibas0=nbask(1)/nrel*(irl-1)
	    call zgemm('n','n',nbask(iq0)/nrel,n,nbask(1)/nrel,
     &	             (1.d0,0.d0),s,nbndtm,ev_bnd(ibas0+1,1,1,ispin),
     &               nbasmpw,(0.d0,0.d0),tmp,nbndtm)
          ibas0=nbask(iq0)/nrel*(irl-1)
	    call zgemm('c','n',n1,n,nbask(iq0)/nrel,(1.d0,0.d0),
     &	             anew(ibas0+1,1),nbasmpw,tmp,nbndtm,(1.d0,0.d0),
     &               bbm(1,1,ipb),n1)
        enddo
      enddo   !!! over id
	deallocate(anew,znew,s,tmp)
      end