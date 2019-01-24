	subroutine sym_vmi_left(ig,f0,f)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
	implicit none
	integer, intent(in) :: ig
	complex*16, intent(in) :: f0(n_pbmt)
	complex*16, intent(out) :: f(n_pbmt)
	integer :: iatom,jatom,isort,ind0,jnd0,l,ie,m,lm,km,lm0
	real*8 :: y0(2*maxpb+1),y1(2*maxpb+1),t0(2*maxpb+1)
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
	        t0(lm0)=dreal(f0(jnd0+km))
	      enddo
	      call rotate1(t0,y0,l,u(1,ig),1)
c ------------- Imaginary components ----------------------------------
            do m=-l,l
              lm=l*(l+1)+m+1
              km=ind_prod(ie,lm,isort)
	        lm0=l+m+1
	        t0(lm0)=dimag(f0(jnd0+km))
	      enddo
	      call rotate1(t0,y1,l,u(1,ig),1)
            do m=-l,l
              lm=l*(l+1)+m+1
              km=ind_prod(ie,lm,isort)
	        lm0=l+m+1
	        f(ind0+km)=dcmplx(y0(lm0),y1(lm0))
	      enddo
	    enddo  !! over ie
	  enddo   !! over l
	enddo   !! over iatom
	end