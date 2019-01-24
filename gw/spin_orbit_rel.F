      subroutine spin_orbit_rel
      use atom_mod
      use etot_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: l,i,mj,i1,ii,ii1,lmj,isort,iatom,mt,n,ist,li1,jj1,
     &	         li,jj,je,in1,jn1,le1,ie,in,jn,le,mt1,irad,limj,limj1,
     &           km,km1
	real*8 :: betta,sqpi4,dqdall,a,b
	real*8, allocatable :: v(:),dv(:)
      complex*16, allocatable :: gloc(:,:,:,:)
	if(irel/=2) return
	allocate(gloc(maxel,maxel,natom,nspin_0))
	gloc=g_loc_0
	sqpi4=sqrt(4.d0*pi)
	if(ubi=='dft') then
	  allocate(v(0:maxnrad))
	  allocate(dv(0:maxnrad))
	endif
	do isort=1,nsort
	  iatom=iat_1(isort)
        if(ubi=='dft') then
          betta=smt(isort)/(exp(h(isort)*nrad(isort))-1.d0)
          mt=indmt(1,isort,1)
	    v(0:nrad(isort))=v_mt_h(mt:mt+nrad(isort))
     &                    +v_mt_xc(mt:mt+nrad(isort))
	    v(0:nrad(isort))=v(0:nrad(isort))/sqpi4
	    dv=0.d0
	    do i=1,nrad(isort)-1
	      dv(i)=(v(i+1)-v(i-1))/2.d0/h(isort)/(r(i,isort)+betta)
     & 	         +2.d0*z(isort)/r(i,isort)**2
          enddo
          n=nrad(isort)
	    dv(n)=(v(n)-v(n-1))/h(isort)/(r(n,isort)+betta)
     & 	       +2.d0*z(isort)/r(n,isort)**2
          dv=dv/clight**2
        endif
	  do l=0,lmb(isort)
	    ist=-1
          if(l==0) ist=1
	    do ii1=ist,1,2
            li1=(5+ii1)/2+2*(l-1)
            jj1=l+l+ii1
			do ii=ist,1,2
              li=(5+ii)/2+2*(l-1)
              jj=l+l+ii
              do je=1,ntle(l,isort)
			    in1=1
			    if(augm(je,l,isort)/='LOC') in1=2
			    do jn1=1,in1
                  le1=indfun0(jn1,je,li1,isort)
	            do ie=1,ntle(l,isort)
			        in=1
			        if(augm(ie,l,isort)/='LOC') in=2
			        do jn=1,in
                      le=indfun0(jn,ie,li,isort)
                      if(ubi=='dft') then
	                  mt=ind_wf(le,isort)
	                  mt1=ind_wf(le1,isort)
                        do irad=0,nrad(isort)
                          v(irad)=gfun(mt+irad,1)*gfun(mt1+irad,1)
                          v(irad)=v(irad)*dv(irad)*dr(irad,isort)
     &                           *r(irad,isort)
                        enddo
                        so_const(le1,le,isort)=dqdall(h(isort),v,
     &                                                nrad(isort))
                      endif
                      a=0.d0
                      do mj=-jj,jj,2
                        call getlimj(limj,l,ii,mj,li,1)
                        km=indbasa(jn,ie,limj,isort)
                        if(abs(mj)>jj1) cycle
                        call getlimj(limj1,l,ii1,mj,li1,1)
                        km1=indbasa(jn1,je,limj1,isort)
                        lmj=l*(l+2)+(mj+3)/2
                        i=(ii+3)/2
                        i1=(ii1+3)/2
                        if(iabs(mj-1)<=l+l) then
                          a=a+(dfloat(mj)/2-0.5d0)*ckg2(i1,2,lmj)
     &                     *ckg2(i,2,lmj)*real(gloc(km,km1,iatom,1))
                          b=sqrt((l+dfloat(mj)/2+0.5d0)
     &                          *(l-dfloat(mj)/2+0.5d0))
                          a=a+b*ckg2(i1,2,lmj)*ckg2(i,1,lmj)
     &                         *real(gloc(km,km1,iatom,1))
                        endif
                        if(iabs(mj+1)<=l+l) then
                          a=a-(dfloat(mj)/2+0.5d0)*ckg2(i1,1,lmj)
     &                     *ckg2(i,1,lmj)*real(gloc(km,km1,iatom,1))
                          b=sqrt((l+dfloat(mj)/2+0.5d0)
     &                          *(l-dfloat(mj)/2+0.5d0))
                          a=a+b*ckg2(i1,1,lmj)*ckg2(i,2,lmj)
     &                         *real(gloc(km,km1,iatom,1))
                        endif
                      enddo   !! over mj
                      e_spin_orbit=e_spin_orbit+a*nhsort(isort)
     &                                         *so_const(le1,le,isort)
                    enddo  !! over jn
                  enddo  !! over ie
                enddo  !! over jn1
              enddo   !! over je
            enddo   !! over ii
          enddo   !! over ii1
        enddo !! over l
      enddo   !! over isort
	if(ubi=='dft') deallocate(v,dv)
	deallocate(gloc)
      end