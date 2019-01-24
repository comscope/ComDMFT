      subroutine set_h0(b,a,nbas,ndim,yr,yi,ind_y,ar0,ai0,br0,bi0,ar1,
     &                  ai1,br1,bi1,ispin)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: nbas
      real*8, intent(in) :: yr(max_bs,nbas),yi(max_bs,nbas)
      integer :: ind_y(2*limlb,natom),ndim,ispin,
     &	         iatom,isort,llmm,ll,ll1,mm,ie1,jbas,lm,l,l1,
     &           m,ie,ibas,li1,ist1,i1,jj1,mj1,li,jj,mj,jatom,
     &           jsort,n,ind,ie0,iba,nbndt,ist,i,ind1,ii,ispin0,ig0,ig,
     &           iatom0,ndimb,ndimb2,lget,jnd
      real*8 :: ar0(nbas,nlcloc),ai0(nbas,nlcloc),hsym(2,2),
     &          br0(nbas,nlcloc),bi0(nbas,nlcloc),ar1(nbas,nbas),
     &          ai1(nbas,nbas),br1(nbas,nbas),bi1(nbas,nbas),
     &          arr,brr,smtcon,drr,di,bbr,
     &          bbi,aar,aai,ffh
      complex*16 :: a(ndim,ndim),b(ndim,ndim)
      complex*16, allocatable :: yc(:,:),zc(:,:,:),tmp(:,:)
	nbndt=nbas/nrel
	ar0=0.d0
	ai0=0.d0
	br0=0.d0
	bi0=0.d0
c --------- ME LO(lo)-LO(lo) (spherical) between local functions -------
      do iatom=1,natom
        isort=is(iatom)
        if(irel.ne.2) then
	    ispin0=ispin
	    iatom0=iat_1(isort)
	    do ig0=1,ngroup
	      if(iatom0==ip(iatom,ig0)) then
	        ig=ig0
	        exit
	      endif
	    enddo
	    llmm=0
          do ll=0,lmb(isort)
            ll1=ll+1
            do mm=-ll,ll
	        llmm=llmm+1
              do ie1=1,ntle(ll,isort)
                if(augm(ie1,ll,isort)=='LOC') then
	            jbas=indbas(ie1,llmm,iatom)  !! the position in basis
	            lm=0
                  do l=0,lmb(isort)
                    l1=l+1
                    do m=-l,l
	                lm=lm+1
                      do ie=1,ntle(l,isort)
                        if(augm(ie,l,isort)=='LOC') then
                          ibas=indbas(ie,lm,iatom)
                          if(ll.eq.l.and.mm.eq.m) then
                            arr=0.5d0*(ffhmt(1,1,ie,ie1,l1,isort,ispin0)
     &				                 +ffhmt(1,1,ie1,ie,l1,isort,ispin0))
                            brr=ffsmt(1,1,ie,ie1,l1,isort,ispin0)
                            a(ibas,jbas)=dcmplx(arr,0.d0)
                            b(ibas,jbas)=dcmplx(brr,0.d0)
                          endif
                        endif
                      enddo
                    enddo
                  enddo
                endif
              enddo   !! over ie1
            enddo  !! over mm
          enddo
        else if(irel.eq.2) then
          li1=0
	    llmm=0
          do ll=0,lmb(isort)
            ist1=-1
            if(ll.eq.0) ist1=1
            do i1=ist1,1,2
              li1=li1+1
              jj1=ll+ll+i1
              do mj1=-jj1,jj1,2
	          llmm=llmm+1
                do ie1=1,ntle(ll,isort)
                  if(augm(ie1,ll,isort)=='LOC') then
                    jbas=indbas(ie1,llmm,iatom)
                    li=0
	              lm=0
                    do l=0,lmb(isort)
                      ist=-1
                      if(l.eq.0) ist=1
                      do i=ist,1,2
                        li=li+1
                        jj=l+l+i
                        do mj=-jj,jj,2
	                    lm=lm+1
                          do ie=1,ntle(l,isort)
                            if(augm(ie,l,isort)=='LOC') then
                              ibas=indbas(ie,lm,iatom)
                              if(li1.eq.li.and.mj.eq.mj1) then
                                arr=0.5d0*
     &                           (ffhmt(1,1,ie,ie1,li,isort,ispin)
     &  	 	                 +ffhmt(1,1,ie1,ie,li,isort,ispin))
                                brr=ffsmt(1,1,ie,ie1,li,isort,ispin)
                                a(ibas,jbas)=dcmplx(arr,0.d0)
                                b(ibas,jbas)=dcmplx(brr,0.d0)
                              endif
                            endif
                          enddo
                        enddo
                      enddo
                    enddo
                  endif
                enddo
              enddo
            enddo
          enddo
        endif
      enddo   !!! over iatom
c ----------- PW-LO(lo) spherical matrix elements ------------------------
      do jatom=1,natom
        jsort=is(jatom)
	  smtcon=0.5d0*smt(jsort)**2
	  if(irel.ne.2) then
 	    ispin0=ispin
	    iatom0=iat_1(jsort)
	    do ig0=1,ngroup
	      if(iatom0==ip(jatom,ig0)) then
	        ig=ig0
	        exit
	      endif
	    enddo
          n=(lmb(jsort)+1)**2
          do l=0,lmb(jsort)
            ie=n_e_apw(l,jsort)
            if(ie==0) cycle
            l1=l+1
            do m=-l,l
              lm=l*l1+m+1
              ind=ind_y(lm,jatom)
              ind1=ind_y(n+lm,jatom)
              do ie0=1,ntle(l,jsort)  !! over local orbitals
                if(augm(ie0,l,jsort)=='LOC') then
                  jbas=indbas(ie0,lm,jatom)
                  do ibas=nlcloc+1,ndim
                    iba=ibas-nlcloc
	              drr=yr(ind,iba)
	              di=-yi(ind,iba)
                    bbr=drr*ffsmt(1,1,ie,ie0,l1,jsort,ispin0)
                    bbi=di*ffsmt(1,1,ie,ie0,l1,jsort,ispin0)
                    aar=0.5d0*drr*(ffhmt(1,1,ie,ie0,l1,jsort,ispin0)
     &			                  +ffhmt(1,1,ie0,ie,l1,jsort,ispin0))
                    aai=0.5d0*di*(ffhmt(1,1,ie,ie0,l1,jsort,ispin0)
     &			 	             +ffhmt(1,1,ie0,ie,l1,jsort,ispin0))
c --------- Standard LAPW calculation for slope-typefunction --------------------
	              drr=yr(ind1,iba)
	              di=-yi(ind1,iba)
                    ffh=0.5d0*(ffhmt(2,1,ie,ie0,l1,jsort,ispin0)
     &                        +ffhmt(1,2,ie0,ie,l1,jsort,ispin0))
                    aar=aar+drr*ffh
                    bbr=bbr+drr*ffsmt(2,1,ie,ie0,l1,jsort,ispin0)
                    aai=aai+di*ffh
                    bbi=bbi+di*ffsmt(2,1,ie,ie0,l1,jsort,ispin0)
                    ar0(iba,jbas)=aar
                    br0(iba,jbas)=bbr
                    ai0(iba,jbas)=aai
                    bi0(iba,jbas)=bbi
                  enddo
                endif
              enddo   !!! over ie0
            enddo   !!! over m
          enddo   !!! over l
        else if(irel.eq.2) then
          n=2*(lmb(jsort)+1)**2
          li=0
          do l=0,lmb(jsort)
            ie=n_e_apw(l,jsort)
            if(ie==0) cycle
            ist=-1
            if(l.eq.0) ist=1
            do i=ist,1,2
              li=li+1
              jj=l+l+i
              do mj=-jj,jj,2
                call getlimj(lm,l,i,mj,li,1)
                ind=ind_y(lm,jatom)
                ind1=ind_y(n+lm,jatom)
                do ie0=1,ntle(l,jsort)  !! over local orbitals
                  if(augm(ie0,l,jsort)=='LOC') then
                    jbas=indbas(ie0,lm,jatom)
                    do ibas=nlcloc+1,ndim
                      iba=ibas-nlcloc
	                drr=yr(ind,iba)
	                di=-yi(ind,iba)
                      bbr=drr*ffsmt(1,1,ie,ie0,li,jsort,ispin)
                      bbi=di*ffsmt(1,1,ie,ie0,li,jsort,ispin)
                      aar=0.5d0*drr*(ffhmt(1,1,ie,ie0,li,jsort,ispin)
     &			                    +ffhmt(1,1,ie0,ie,li,jsort,ispin))
                      aai=0.5d0*di*(ffhmt(1,1,ie,ie0,li,jsort,ispin)
     &			                   +ffhmt(1,1,ie0,ie,li,jsort,ispin))
c --------- Standard LAPW calculation for slope-typefunction --------------------
	                drr=yr(ind1,iba)
	                di=-yi(ind1,iba)
                      ffh=0.5d0*(ffhmt(2,1,ie,ie0,li,jsort,ispin)
     &                          +ffhmt(1,2,ie0,ie,li,jsort,ispin))
                      aar=aar+drr*ffh
                      bbr=bbr+drr*ffsmt(2,1,ie,ie0,li,jsort,ispin)
                      aai=aai+di*ffh
                      bbi=bbi+di*ffsmt(2,1,ie,ie0,li,jsort,ispin)
                      ar0(iba,jbas)=aar
                      br0(iba,jbas)=bbr
                      ai0(iba,jbas)=aai
                      bi0(iba,jbas)=bbi
                    enddo
                  endif
                enddo   !!! over ie0
              enddo   !!! over mj
            enddo   !!! over i
          enddo   !!! over l
        endif
      enddo   !!! over jatom
c ----------- PW-PW spherical+interstitial matrix elements -----------
      allocate(tmp(nbas,nbas))
      do iatom=1,natom
        isort=is(iatom)
        ndimb=nrel*(lmb(isort)+1)**2
        ind=ind_y(1,iatom)
        ndimb2=2*ndimb
        allocate(yc(ndimb2,nbas))
        allocate(zc(ndimb2,nbas,2))
        do i=1,nbas
          yc(:,i)=dcmplx(yr(ind:ind+ndimb2-1,i),yi(ind:ind+ndimb2-1,i))
        enddo
        zc=(0.d0,0.d0)
        do lm=1,ndimb
	    if(irel.le.1) then
	      l=lget(lm)
	      li=l+1
	    else
	      call getlimj(lm,l,i,m,li,0)
	    endif
          ie=n_e_apw(l,isort)
          hsym=ffhmt(1:2,1:2,ie,ie,li,isort,ispin)
          hsym(1,2)=0.5d0*(hsym(1,2)+hsym(2,1))
          hsym(2,1)=hsym(1,2)
          do jj=1,2
            jnd=(jj-1)*ndimb+lm
            do ii=1,2
              ind=(ii-1)*ndimb+lm
              zc(ind,:,1)=zc(ind,:,1)+ffsmt(ii,jj,ie,ie,li,isort,ispin)
     &                               *yc(jnd,:)
              zc(ind,:,2)=zc(ind,:,2)+hsym(ii,jj)*yc(jnd,:)
            enddo
          enddo
        enddo
        call zgemm('c','n',nbas,nbas,ndimb2,(1.d0,0.d0),yc,ndimb2,
     &             zc(1,1,1),ndimb2,(0.d0,0.d0),tmp,nbas)
        br1=br1+real(tmp)
        bi1=bi1+imag(tmp)
        call zgemm('c','n',nbas,nbas,ndimb2,(1.d0,0.d0),yc,ndimb2,
     &             zc(1,1,2),ndimb2,(0.d0,0.d0),tmp,nbas)
        ar1=ar1+real(tmp)
        ai1=ai1+imag(tmp)
        deallocate(yc,zc)
      enddo    !! over iatom
      deallocate(tmp)
      end