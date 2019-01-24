      subroutine split_magn(key)
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
      logical :: lb
      integer :: isort,it,iatom,lm,l,lget,li,je,in1,jn1,kmj,ie,in,jn,km,
     &           l1,ist1,ist,ii1,li1,le1,mt1,jj1,ii,le,mt,jj,ir,mj,
     &           lm1,mj1
      real*8 :: s,c2,ss,dqdall,uu,b(3),uu1
      real*8, allocatable :: work(:)
      if(ubi=='dft'.and.key/=2) return
      if(ubi/='dft'.and.key/=0) return
      lb=.false.
      if(ubi=='dft') then
        it=iter_h_ext/10000
        if(iter<=it) lb=.true.
      else if(ubi==' hf') then
        it=mod(iter_h_ext,10000)/1000
        if(iter-iter_dft<=it) lb=.true.
      else if(ubi==' gw') then
        it=mod(iter_h_ext,1000)/100
        if(iter-iter_dft-iter_hf<=it) lb=.true.
      else if(ubi==' qp') then
        it=mod(iter_h_ext,100)/10
        if(iter-iter_dft-iter_hf-iter_gw<=it) lb=.true.
      else if(ubi=='gwg') then
        it=mod(iter_h_ext,10)
        if(iter-iter_dft-iter_hf-iter_gw-iter_qp<=it) lb=.true.
      endif
	if(.not.lb) return
	c2=clight**2
	if(irel==2) allocate(work(0:maxnrad))
	do iatom=1,natom
	  isort=is(iatom)
	  b=-b_ext*magn_shift(isort)
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
	              s=magn_shift(isort)
     &               *ffsmt(jn,jn1,ie,je,li,isort,1)
                    pv(km,kmj,iatom,1)=pv(km,kmj,iatom,1)+s
	              s=magn_shift(isort)
     &               *ffsmt(jn,jn1,ie,je,li,isort,nspin_1)
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
c ----- Big components contribution ---------------------------
                          if(l==l1) then
                            do ir=0,nrad(isort)
                              work(ir)=gfun(mt+ir,1)*gfun(mt1+ir,1)
     &                                *dr(ir,isort)*r(ir,isort)**2
                            enddo
                            s=dqdall(h(isort),work,nrad(isort))
c   -------  B(z) contribution -------------------------------                            
                            do mj=-jj,jj,2
                              if(mj<-jj1.or.mj>jj1) cycle
		                      call getlimj(lm,l,ii,mj,li,1)
		                      call getlimj(lm1,l1,ii1,mj,li1,1)
	                        km=indbasa(jn,ie,lm,isort)
	                        kmj=indbasa(jn1,je,lm1,isort)
                              uu=dfloat(mj)/dfloat(l+l+1)
                              ss=0.d0
                              if(ii==-1.and.ii1==-1) then
                                uu1=0.5d0*(1.d0+uu)
                                uu=0.5d0*(1.d0-uu)
                                if(mj+1<=l+l) ss=ss-uu1
                                if(mj-1>=-l-l) ss=ss+uu
                              else if(ii/=ii1) then
                                uu=0.5d0*sqrt(1.d0-uu*uu)
                                if(mj+1<=l+l) ss=ss-uu
                                if(mj-1>=-l-l) ss=ss-uu
                              else if(ii==1.and.ii1==1) then
                                uu1=0.5d0*(1.d0+uu)
                                uu=0.5d0*(1.d0-uu)
                                if(mj+1<=l+l) ss=ss-uu
                                if(mj-1>=-l-l) ss=ss+uu1
                              endif
                              pv(km,kmj,iatom,1)=pv(km,kmj,iatom,1)
     &                                          +s*ss*b(3)
                            enddo
c   -------  B(-) contribution -------------------------------
                            do mj=-jj,jj,2
                              if(mj-1<-l-l) cycle
                              if(mj-1>l+l) cycle
                              mj1=mj-2
                              if(mj1+1>l+l) cycle
                              if(mj1<-jj1) cycle
                              if(mj1>jj1) cycle
		                      call getlimj(lm,l,ii,mj,li,1)
		                      call getlimj(lm1,l1,ii1,mj1,li1,1)
	                        km=indbasa(jn,ie,lm,isort)
	                        kmj=indbasa(jn1,je,lm1,isort)
                              uu=dfloat(mj)/dfloat(l+l+1)
                              uu1=dfloat(mj1)/dfloat(l+l+1)
                              ss=0.0d0
                              if(ii==-1.and.ii1==-1) then
                                ss=-0.5d0*sqrt((1.d0-uu)*(1.d0+uu1))
                              else if(ii==-1.and.ii1==1) then
                                ss=-0.5d0*sqrt((1.d0-uu)*(1.d0-uu1))
                              else if(ii==1.and.ii1==-1) then
                                ss=0.5d0*sqrt((1.d0+uu)*(1.d0+uu1))
                              else if(ii==1.and.ii1==1) then
                                ss=0.5d0*sqrt((1.d0+uu)*(1.d0-uu1))
                              endif
                              pv(km,kmj,iatom,1)=pv(km,kmj,iatom,1)
     &                                          +s*ss*b(1)
                              pvj(km,kmj,iatom)=pvj(km,kmj,iatom)
     &                                         -s*ss*b(2)
                            enddo
c   -------  B(+) contribution -------------------------------                            
                            do mj=-jj,jj,2
                              if(mj+1>l+l) cycle
                              if(mj+1<-l-l) cycle
                              mj1=mj+2
                              if(mj1-1<-l-l) cycle
                              if(mj1>jj1) cycle
                              if(mj1<-jj1) cycle
		                      call getlimj(lm,l,ii,mj,li,1)
		                      call getlimj(lm1,l1,ii1,mj1,li1,1)
	                        km=indbasa(jn,ie,lm,isort)
	                        kmj=indbasa(jn1,je,lm1,isort)
                              uu=dfloat(mj)/dfloat(l+l+1)
                              uu1=dfloat(mj1)/dfloat(l+l+1)
                              ss=0.0d0
                              if(ii==-1.and.ii1==-1) then
                                ss=-0.5d0*sqrt((1.d0+uu)*(1.d0-uu1))
                              else if(ii==-1.and.ii1==1) then
                                ss=0.5d0*sqrt((1.d0+uu)*(1.d0+uu1))
                              else if(ii==1.and.ii1==-1) then
                                ss=-0.5d0*sqrt((1.d0-uu)*(1.d0-uu1))
                              else if(ii==1.and.ii1==1) then
                                ss=0.5d0*sqrt((1.d0-uu)*(1.d0+uu1))
                              endif
                              pv(km,kmj,iatom,1)=pv(km,kmj,iatom,1)
     &                                          +s*ss*b(1)
                              pvj(km,kmj,iatom)=pvj(km,kmj,iatom)
     &                                         +s*ss*b(2)
                            enddo
                          endif
c ----- Small components contribution ---------------------------
                          if(l+ii==l1+ii1+1000) then
                            do ir=0,nrad(isort)
                              work(ir)=gfund(mt+ir,1)*gfund(mt1+ir,1)
     &                                *dr(ir,isort)*r(ir,isort)**2
                            enddo
                            s=magn_shift(isort)
     *                       *dqdall(h(isort),work,nrad(isort))/c2
                            do mj=-jj,jj,2
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
	call check_pv_sym
      end
