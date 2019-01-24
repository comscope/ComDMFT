      subroutine sym_me
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
	integer :: isort,iatom,ja,ig,jatom,j,l,lp1,ie,m,lm,km,lm0,li,i,
     &	         ist,ii,jj,mj,ispin,jspin,in,jn
      real*8, allocatable :: y0(:),t0(:),pa1(:,:)
	complex*16, allocatable :: y0j(:),pa1j(:,:),t0j(:)
	allocate(y0(2*maxb+1),t0(2*maxb+1),pa1(maxel,maxel))
	allocate(y0j(2*maxb+2),pa1j(maxel,maxel),t0j(2*maxb+2))
c ------ start loop over atoms in class ----------
      do ispin=1,nspin
	  do isort=1,nsort
	    iatom=iat_1(isort)
          do ja=iatom+1,natom     !!  over atoms of the given class
            do ig=1,ngroup
              jatom=ip(ja,ig)
	        jspin=ispin
	        if(jatom.eq.iatom) then
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
     &					  	                  pvj(km,j,iatom))
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
	                      pv(i,km,ja,jspin)=y0(lm0)
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
	                        pv(i,km,ja,jspin)=dreal(y0j(lm0))
	                        pvj(i,km,ja)=dimag(y0j(lm0))
	                      enddo
	                    enddo
	                  enddo   !! over ie
	                enddo   !! over ii
	              enddo   !! over l
	            endif
	          enddo  !! over i
	          goto 2
	        endif
	      enddo   !! over ig
2           continue
          enddo   !! over ja
        enddo  !! over isort
	enddo   !! over ispin
	deallocate(y0,t0,pa1)
	deallocate(y0j,pa1j,t0j)
      end