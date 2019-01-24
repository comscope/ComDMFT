      subroutine check_pv_sym
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer :: isort,iatom,ig,j,l,lp1,ie,m,lm,km,lm0,li,i,
     &	         ist,ii,jj,mj,ispin,in,jn,ngat
	real*8 :: diff,d
      real*8, allocatable :: y0(:),t0(:),pa1(:,:)
	complex*16, allocatable :: y0j(:),pa1j(:,:),t0j(:),tmp(:,:)
	allocate(y0(2*maxb+1),t0(2*maxb+1),pa1(maxel,maxel))
	allocate(y0j(2*maxb+2),pa1j(maxel,maxel),t0j(2*maxb+2))
	allocate(tmp(maxel,maxel))
      do ispin=1,nspin
	  do iatom=1,natom
	    isort=is(iatom)
	    tmp=(0.d0,0.d0)
	    ngat=0
          do ig=1,ngroup
            if(ip(iatom,ig)/=iatom) cycle
            ngat=ngat+1
	      do j=1,lfunm(isort)
              if(irel<=1) then
	          do l=0,lmb(isort)
                  lp1=l+1
                  do ie=1,ntle(l,isort)
			        in=1
			        if(augm(ie,l,isort)/='LOC') in=2
			        do jn=1,in
                      do m=-l,l
                        lm=l*(l+1)+m+1
                        km=indbasa(jn,ie,lm,isort)
	                  lm0=l+m+1
	                  t0(lm0)=pv(km,j,iatom,ispin)
	                enddo
	                call rotate1(t0,y0,l,u(1,ig),1)
                      do m=-l,l
                        lm=l*(l+1)+m+1
                        km=indbasa(jn,ie,lm,isort)
	                  lm0=l+m+1
	                  pa1(km,j)=y0(lm0)
	                enddo
	              enddo
	            enddo
	          enddo   !! over l
	        else if(irel.eq.2) then
                li=0
                do l=0,lmb(isort)
                  ist=-1
                  if(l.eq.0) ist=1
                  do ii=ist,1,2
                    li=li+1
                    jj=l+l+ii
                    do ie=1,ntle(l,isort)
			          in=1
			          if(augm(ie,l,isort)/='LOC') in=2
			          do jn=1,in
	                  lm0=0
                        do mj=-jj,jj,2
                          lm0=lm0+1
			              call getlimj(lm,l,ii,mj,li,1)
                          km=indbasa(jn,ie,lm,isort)
	                    t0j(lm0)=dcmplx(pv(km,j,iatom,ispin),
     &				  	                  pvj(km,j,iatom))
	                  enddo
	                  call rotate2(t0j,y0j,l,ii,uj(1,ig),1,1)
                        do mj=-jj,jj,2
			              call getlimj(lm,l,ii,mj,li,1)
                          km=indbasa(jn,ie,lm,isort)
	                    lm0=(jj+mj)/2+1
	                    pa1j(km,j)=y0j(lm0)
	                  enddo
	                enddo
	              enddo
	            enddo
	          enddo
	        endif
	      enddo
	      do i=1,lfunm(isort)
              if(irel<=1) then
		        do l=0,lmb(isort)
                  lp1=l+1
                  do ie=1,ntle(l,isort)
			        in=1
			        if(augm(ie,l,isort)/='LOC') in=2
			        do jn=1,in
                      do m=-l,l
                        lm=l*(l+1)+m+1
                        km=indbasa(jn,ie,lm,isort)
	                  lm0=l+m+1
	                  t0(lm0)=pa1(i,km)
	                enddo
	                call rotate1(t0,y0,l,u(1,ig),1)
                      do m=-l,l
                        lm=l*(l+1)+m+1
                        km=indbasa(jn,ie,lm,isort)
	                  lm0=l+m+1
	                  tmp(i,km)=tmp(i,km)+y0(lm0)
	                enddo
	              enddo
	            enddo
	          enddo
	        else if(irel.eq.2) then
                li=0
                do l=0,lmb(isort)
                  ist=-1
                  if(l.eq.0) ist=1
                  do ii=ist,1,2
                    li=li+1
                    jj=l+l+ii
                    do ie=1,ntle(l,isort)
			          in=1
			          if(augm(ie,l,isort)/='LOC') in=2
			          do jn=1,in
	                  lm0=0
                        do mj=-jj,jj,2
                          lm0=lm0+1
			              call getlimj(lm,l,ii,mj,li,1)
                          km=indbasa(jn,ie,lm,isort)
	                    t0j(lm0)=pa1j(i,km)
	                  enddo
	                  call rotate2(t0j,y0j,l,ii,uj(1,ig),1,2)
                        do mj=-jj,jj,2
			              call getlimj(lm,l,ii,mj,li,1)
                          km=indbasa(jn,ie,lm,isort)
	                    lm0=(jj+mj)/2+1
	                    tmp(i,km)=tmp(i,km)+y0j(lm0)
	                  enddo
	                enddo
	              enddo   !! over ie
	            enddo   !! over ii
	          enddo   !! over l
	        endif
	      enddo  !! over i
	    enddo   !! over ig
	    if(ngat/=0) then
	      tmp=tmp/ngat
	      diff=0.d0
	      do j=1,lfunm(isort)
	        do i=1,lfunm(isort)
	          if(irel==2) d=abs(tmp(i,j)
     &	                       -dcmplx(pv(i,j,iatom,1),pvj(i,j,iatom)))
	          if(irel/=2) d=abs(tmp(i,j)
     &	                       -dcmplx(pv(i,j,iatom,ispin),0.d0))
c                if(d>1.d-4) print *,i,j,d
                diff=max(diff,d)
              enddo
            enddo
            if(maswrk) write(iun,*)' CHECK_PV: ',iatom,diff
          endif
        enddo  !! over iatom
	enddo   !! over ispin
	deallocate(y0,t0,pa1,tmp,y0j,pa1j,t0j)
      end