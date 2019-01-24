	subroutine sym_z_0(znew,k,z_bnd0,neig,ig,pt)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
	implicit none
	integer, intent(in) :: k,neig,ig
	real*8, intent(in) :: pt(3)
	complex*16, intent(in) :: z_bnd0(nfun,neig)
	complex*16,intent(out) :: znew(nfun,neig)
	integer :: ibnd,iatom,jatom,isort,ind0,jnd0,l,ie,m,k0,
     &	         lm,km,lm0,li,ist,jj,ii,mj,ind1,in,jn
	real*8 :: pi2,fas,v(3)
	complex*16 :: cf
	real*8, allocatable :: y0(:),y1(:),t0(:)
	complex*16, allocatable :: t0j(:),y0j(:)
	k0=i_kref(k)
	if(k==k0) then
	  znew=z_bnd0
	  return
	endif
	if(irel<=1) allocate(y0(2*maxb+1),y1(2*maxb+1),t0(2*maxb+1))
	if(irel==2) allocate(t0j(2*maxb+2),y0j(2*maxb+2))
	pi2=pi+pi
	do ibnd=1,neig
	  do iatom=1,natom
	    jatom=ip(iatom,ig)
	    isort=is(iatom)
	    ind0=io_lem(iatom)-1   !! zero index for atom iatom
	    jnd0=io_lem(jatom)-1
          if(irel<=1) then
	      do l=0,lmb(isort)
              do ie=1,ntle(l,isort)
	          in=1
		        if(augm(ie,l,isort)/='LOC') in=2
		        do jn=1,in
c ------------- Real components ----------------------------------
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    km=indbasa(jn,ie,lm,isort)
	              lm0=l+m+1
	              t0(lm0)=dreal(z_bnd0(jnd0+km,ibnd))
	            enddo
	            call rotate1(t0,y0,l,u(1,ig),1)
c ------------- Imaginary components ----------------------------------
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    km=indbasa(jn,ie,lm,isort)
	              lm0=l+m+1
	              t0(lm0)=dimag(z_bnd0(jnd0+km,ibnd))
	            enddo
	            call rotate1(t0,y1,l,u(1,ig),1)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    km=indbasa(jn,ie,lm,isort)
	              lm0=l+m+1
	              znew(ind0+km,ibnd)=dcmplx(y0(lm0),y1(lm0))
	            enddo
	          enddo
	        enddo  !! over ie
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
		    	      t0j(lm0)=z_bnd0(jnd0+km,ibnd)
	              enddo
	              call rotate2(t0j,y0j,l,ii,uj(1,ig),1,1)
                    do mj=-jj,jj,2
			          call getlimj(lm,l,ii,mj,li,1)
                      km=indbasa(jn,ie,lm,isort)
	                lm0=(jj+mj)/2+1
	                znew(ind0+km,ibnd)=y0j(lm0)
	              enddo
	            enddo
	          enddo
	        enddo
	      enddo
	    endif
c ---------------  VII-30-21 -----------------------------------------
          ind1=0
          if(iatom<natom) ind1=io_lem(iatom+1)-1 !! last index for iatom
	    if(iatom==natom) ind1=nfun !! last index for atom iatom
	    v=tshift(:,iatom,ig)
 	    fas=pi2*dot_product(pt,v)
	    cf=dcmplx(cos(fas),sin(fas))
	    znew(ind0+1:ind1,ibnd)=cf*znew(ind0+1:ind1,ibnd)
	  enddo   !! over iatom
	enddo  !! over ibnd
	if(irel<=1) deallocate(y0,y1,t0)
	if(irel==2) deallocate(t0j,y0j)
	end
