	subroutine sym_z(znew,k,z_bnd0,neig)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
	implicit none
	integer, intent(in) :: k,neig
	complex*16, intent(in) :: z_bnd0(nfun,neig)
	complex*16,intent(out) :: znew(nfun,neig)
	integer :: ibnd,iatom,jatom,isort,ind0,jnd0,l,ie,m,k0,lm,km,lm0,
     &	         in,jn,ns,ig
	real*8 :: pi2,phase
	complex*16 :: cc
	real*8, allocatable :: y0(:),y1(:),t0(:)
	pi2=pi+pi
	allocate(y0(2*maxb+1),y1(2*maxb+1),t0(2*maxb+1))
	k0=i_kref(k)
	if(k==k0) then
	  znew=z_bnd0
	  return
	endif
	ns=0
	znew=(0.d0,0.d0)
	do ig=1,ngroup
	  if(ip_k(k,ig)/=k0) cycle
	  ns=ns+1
	  do iatom=1,natom
	    jatom=ip(iatom,ig)
	    isort=is(iatom)
	    ind0=io_lem(iatom)-1   !! zero index for atom iatom
	    jnd0=io_lem(jatom)-1
 	    phase=pi2*dot_product(pnt(:,k),tshift(:,iatom,ig))
	    cc=dcmplx(cos(phase),sin(phase))
	    do l=0,lmb(isort)
            do ie=1,ntle(l,isort)
	        in=1
		      if(augm(ie,l,isort)/='LOC') in=2
	        do ibnd=1,neig
		        do jn=1,in
c ----------- Real components ----------------------------------
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
	              znew(ind0+km,ibnd)=znew(ind0+km,ibnd)
     &	                              +cc*dcmplx(y0(lm0),y1(lm0))
	            enddo
	          enddo
	        enddo   !! over ibnd
	      enddo  !! over ie
	    enddo   !! over l
	  enddo   !! over iatom
	enddo  !! over ig
	deallocate(y0,y1,t0)
	znew=znew/ns
	end