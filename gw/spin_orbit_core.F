      subroutine spin_orbit_core
      use atom_mod
      use etot_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: l,i,mj,i1,ii,ii1,lmj,isort,iatom,mt,n,li1,jj1,li,jj,je,
     &	         ie,irad,limj,limj1,jcore,l1,konf1,icore,konf
	real*8 :: betta,sqpi4,dqdall,a,b,so
	real*8, allocatable :: v(:),dv(:)
	if(irel_core/=2) return
	sqpi4=sqrt(4.d0*pi)
	allocate(v(0:maxnrad))
	allocate(dv(0:maxnrad))
	do isort=1,nsort
	  iatom=iat_1(isort)
        betta=smt(isort)/(exp(h(isort)*nrad(isort))-1.d0)
        mt=indmt(1,isort,1)
	  v(0:nrad(isort))=v_mt_h(mt:mt+nrad(isort))
c     &                  +v_mt_xc(mt:mt+nrad(isort))
	  v(0:nrad(isort))=v(0:nrad(isort))/sqpi4
	  dv=0.d0
	  do i=1,nrad(isort)-1
	    dv(i)=(v(i+1)-v(i-1))/2.d0/h(isort)/(r(i,isort)+betta)
     & 	       +2.d0*z(isort)/r(i,isort)**2
        enddo
        n=nrad(isort)
	  dv(n)=(v(n)-v(n-1))/h(isort)/(r(n,isort)+betta)
     & 	     +2.d0*z(isort)/r(n,isort)**2
        dv=dv/clight**2	
	  jcore=0
        do li1=1,2*lmb(isort)+1
	    l1=li1/2
	    ii1=li1-2*l1
	    if(ii1.eq.0) ii1=-1
          do konf1=l1+1,konfig(l1,isort)-1
            jcore=jcore+1
	      je=indcor(jcore,1,isort)
	      icore=0
            do li=1,2*lmb(isort)+1
	        l=li/2
	        ii=li-2*l
	        if(ii.eq.0) ii=-1
              do konf=l+1,konfig(l,isort)-1
                icore=icore+1
	          ie=indcor(icore,1,isort)
                if(l1/=l) cycle
                jj1=l+l+ii1
                jj=l+l+ii
                v(0)=0.d0
                do irad=1,nrad(isort)
                  v(irad)=pcor(ie+irad)*pcor(je+irad)
                  v(irad)=v(irad)*dv(irad)*dr(irad,isort)/r(irad,isort)
                enddo
                so=dqdall(h(isort),v,nrad(isort))
                a=0.d0
                do mj=-jj,jj,2
                  call getlimj(limj,l,ii,mj,li,1)
                  if(abs(mj)>jj1) cycle
                  call getlimj(limj1,l,ii1,mj,li1,1)
                  lmj=l*(l+2)+(mj+3)/2
                  i=(ii+3)/2
                  i1=(ii1+3)/2
                  if(iabs(mj-1)<=l+l) then
                    a=a+(dfloat(mj)/2-0.5d0)*ckg2(i1,2,lmj)
     &                                      *ckg2(i,2,lmj)
                    b=sqrt((l+dfloat(mj)/2+0.5d0)
     &                    *(l-dfloat(mj)/2+0.5d0))
                    a=a+b*ckg2(i1,2,lmj)*ckg2(i,1,lmj)
                  endif
                  if(iabs(mj+1)<=l+l) then
                    a=a-(dfloat(mj)/2+0.5d0)*ckg2(i1,1,lmj)
     &                                      *ckg2(i,1,lmj)
                    b=sqrt((l+dfloat(mj)/2+0.5d0)
     &                    *(l-dfloat(mj)/2+0.5d0))
                    a=a+b*ckg2(i1,1,lmj)*ckg2(i,2,lmj)
                  endif
                enddo   !! over mj
                e_spin_orbit_core=e_spin_orbit_core+a*nhsort(isort)*so
              enddo   !! over konf
            enddo   !! over li
          enddo   !! over konf1
        enddo !! over li1
      enddo   !! over isort
	deallocate(v,dv)
      end