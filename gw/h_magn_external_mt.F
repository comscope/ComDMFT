      subroutine h_magn_external_mt(key)
c	key=0 - only Hartree part
c	key=1 - Only XC part
c	key=2 - Hartree + XC
	use atom_mod
	use etot_mod
      use manager_mod
      use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: key
      integer :: isort,iatom,lm,l,lget,li,je,in1,jn1,kmj,ie,in,jn,km,
     &           ind0,l1,ist1,ist,ii1,li1,le1,mt1,jj1,ii,le,mt,jj,jj0,
     &           ir,mj,lm1
      real*8 :: s,c2,ss,dqdall,uu
      real*8, allocatable :: work(:)
      if(ubi/='dft') return
      if(key/=2) return
	c2=clight**2
	if(irel==2) allocate(work(0:maxnrad))
	do iatom=1,natom
	  isort=is(iatom)
	  ind0=io_lem(iatom)-1
	  if(irel<=1) then
	    do lm=1,(lmb(isort)+1)**2
		    l=lget(lm)
		    li=l+1
	      do je=1,ntle(l,isort)
			  in1=1
			  if(augm(je,l,isort)/='LOC') in1=2
			  do jn1=1,in1
	          kmj=indbasa(jn1,je,lm,isort)
	          do ie=1,ntle(l,isort)
			      in=1
			      if(augm(ie,l,isort)/='LOC') in=2
			      do jn=1,in
	              km=indbasa(jn,ie,lm,isort)
	              s=b_extval*ffsmt(jn,jn1,ie,je,li,isort,1)
                    pv(km,kmj,iatom,1)=pv(km,kmj,iatom,1)+s
	              s=b_extval*ffsmt(jn,jn1,ie,je,li,isort,nspin_1)
                    pv(km,kmj,iatom,nspin)=pv(km,kmj,iatom,nspin)-s
	            enddo
	          enddo
	        enddo
	      enddo
	    enddo    !! over lm
	  else if(irel==2) then
	    do l1=0,lmb(isort)
            ist1=-1
            if(l1==0) ist1=1
	      do l=0,lmb(isort)
              ist=-1
              if(l==0) ist=1
	        do je=1,ntle(l1,isort)
			    in1=1
			    if(augm(je,l1,isort)/='LOC') in1=2
			    do jn1=1,in1
	            do ie=1,ntle(l,isort)
			        in=1
			        if(augm(ie,l,isort)/='LOC') in=2
			        do jn=1,in
			          do ii1=ist1,1,2
                        li1=(5+ii1)/2+2*(l1-1)
                        le1=indfun0(jn1,je,li1,isort)
	                  mt1=ind_wf(le1,isort)
                        jj1=l1+l1+ii1
			            do ii=ist,1,2
                          li=(5+ii)/2+2*(l-1)
                          le=indfun0(jn,ie,li,isort)
	                    mt=ind_wf(le,isort)
                          jj=l+l+ii
                          jj0=min(jj,jj1)
c ----- Big components contribution ---------------------------
                          if(l==l1) then
                            do ir=0,nrad(isort)
                              work(ir)=gfun(mt+ir,1)*gfun(mt1+ir,1)
     &                                *dr(ir,isort)*r(ir,isort)**2
                            enddo
                            s=-b_extval
     &                        *dqdall(h(isort),work,nrad(isort))
                            do mj=-jj0,jj0,2
		                      call getlimj(lm,l,ii,mj,li,1)
		                      call getlimj(lm1,l1,ii1,mj,li1,1)
	                        km=indbasa(jn,ie,lm,isort)
	                        kmj=indbasa(jn1,je,lm1,isort)
                              uu=dfloat(mj)/dfloat(l+l+1)
                              ss=0.0d0
                              if(ii==-1.and.ii1==-1) ss=-uu
                              if(ii==-1.and.ii1==1) ss=-sqrt(1.d0-uu*uu)
                              if(ii==1.and.ii1==-1) ss=-sqrt(1.d0-uu*uu)
                              if(ii==1.and.ii1==1) ss=uu
                              pv(km,kmj,iatom,1)=pv(km,kmj,iatom,1)+s*ss
                            enddo
                          endif
c ----- Small components contribution ---------------------------
                          if(l+ii==l1+ii1) then
                            do ir=0,nrad(isort)
                              work(ir)=gfund(mt+ir,1)*gfund(mt1+ir,1)
     &                                *dr(ir,isort)*r(ir,isort)**2
                            enddo
                            s=b_extval
     *                       *dqdall(h(isort),work,nrad(isort))/c2
                            do mj=-jj0,jj0,2
		                      call getlimj(lm,l,ii,mj,li,1)
		                      call getlimj(lm1,l1,ii1,mj,li1,1)
	                        km=indbasa(jn,ie,lm,isort)
	                        kmj=indbasa(jn1,je,lm1,isort)
                              uu=dfloat(mj)/dfloat(2*(l+ii)+1)
                              ss=0.0d0
                              if(ii==-1.and.ii1==-1) ss=uu
                              if(ii==-1.and.ii1==1) ss=-sqrt(1.d0-uu*uu)
                              if(ii==1.and.ii1==-1) ss=-sqrt(1.d0-uu*uu)
                              if(ii==1.and.ii1==1) ss=-uu
                              pv(km,kmj,iatom,1)=pv(km,kmj,iatom,1)+s*ss
                            enddo
                          endif
                        enddo  !! over ii
                      enddo  !! over ii1
                    enddo !! over jn
                  enddo  !! over ie
                enddo  !! over jn1
              enddo  !! over je
            enddo   !! over l
          enddo   !! over l1
	  endif
	enddo  !! iatom
	if(irel==2) deallocate(work)
      end
