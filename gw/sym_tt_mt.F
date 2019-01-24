	subroutine sym_tt_mt(tt,k,tt0,nop,ig,pt)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
	implicit none
	integer, intent(in) :: nop,ig,k
	real*8, intent(in) :: pt(3)
	complex*16, intent(in) :: tt0(n_pbtot,nop)
	complex*16,intent(out) :: tt(n_pbmt,nop)
	integer :: ibnd,iatom,jatom,isort,ind0,jnd0,l,ie,m,lm,km,lm0,k0,
     & 	         ind1
	real*8 :: pi2,fas,v(3)
	complex*16 :: cf
	real*8, allocatable :: y0(:),y1(:),t0(:)
	k0=i_kref(k)
	if(k==k0) then
	  tt=tt0(1:n_pbmt,:)
	  return
	endif
	allocate(y0(2*maxpb+1),y1(2*maxpb+1),t0(2*maxpb+1))
	pi2=pi+pi
	do ibnd=1,nop
	  do iatom=1,natom
	    jatom=ip(iatom,ig)
	    isort=is(iatom)
	    ind0=iopb(iatom)-1   !! zero index for atom iatom
	    jnd0=iopb(jatom)-1
	    do l=0,lmpb(isort)
            do ie=1,ntle_pb(l,isort)
c ------------- Real components ----------------------------------
              do m=-l,l
                lm=l*(l+1)+m+1
                km=ind_prod(ie,lm,isort)
	          lm0=l+m+1
	          t0(lm0)=dreal(tt0(jnd0+km,ibnd))
	        enddo
	        call rotate1(t0,y0,l,u(1,ig),1)
c ------------- Imaginary components ----------------------------------
              do m=-l,l
                lm=l*(l+1)+m+1
                km=ind_prod(ie,lm,isort)
	          lm0=l+m+1
	          t0(lm0)=dimag(tt0(jnd0+km,ibnd))
	        enddo
	        call rotate1(t0,y1,l,u(1,ig),1)
              do m=-l,l
                lm=l*(l+1)+m+1
                km=ind_prod(ie,lm,isort)
	          lm0=l+m+1
	          tt(ind0+km,ibnd)=dcmplx(y0(lm0),y1(lm0))
	        enddo
	      enddo  !! over ie
	    enddo   !! over l
c ---------------  Symmetry-19 -----------------------------------------
          ind1=0
          if(iatom<natom) ind1=iopb(iatom+1)-1 !! last index for iatom
	    if(iatom==natom) ind1=n_pbmt !! last index for atom iatom
	    v=tshift(:,iatom,ig)
 	    fas=pi2*dot_product(pt,v)
	    cf=dcmplx(cos(fas),sin(fas))
	    tt(ind0+1:ind1,ibnd)=cf*tt(ind0+1:ind1,ibnd)
	  enddo   !! over iatom
	enddo  !! over ibnd
	deallocate(y0,y1,t0)
	end
