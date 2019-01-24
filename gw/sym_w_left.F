      subroutine sym_w_left(ip0,u0,s,key)
c	key=1 Y(m)=Sum_m1 D_m_m1 X_m1
c	key=2 Y(m)=Sum_m1 D_m1_m X_m1
	use atom_mod
	use solid_mod
      implicit none
	integer :: ip0(natom),iatom,ind0,isort,jatom,jnd0,l,ie,m,lm,km,
     &	         lm0,key
      real*8 :: u0(maxwig)
	complex*16 :: s(n_pbmt)
	real*8, allocatable :: y0(:),y1(:),t0(:)
	complex*16, allocatable :: s0(:)
	allocate(y0(2*maxpb+1),t0(2*maxpb+1),y1(2*maxpb+1),s0(n_pbmt))
	s0=s
      do iatom=1,natom     !!  over atoms
	  ind0=iopb(iatom)-1
	  isort=is(iatom)
        jatom=ip0(iatom)
	  jnd0=iopb(jatom)-1
	  do l=0,lmpb(isort)
          do ie=1,ntle_pb(l,isort)
c ------------- Real components ----------------------------------
            do m=-l,l
              lm=l*(l+1)+m+1
              km=ind_prod(ie,lm,isort)
	        lm0=l+m+1
	        t0(lm0)=dreal(s0(jnd0+km))
	      enddo
	      call rotate1(t0,y0,l,u0,key)
c ------------- Imaginary components ----------------------------------
            do m=-l,l
              lm=l*(l+1)+m+1
              km=ind_prod(ie,lm,isort)
	        lm0=l+m+1
	        t0(lm0)=dimag(s0(jnd0+km))
	      enddo
	      call rotate1(t0,y1,l,u0,key)
c -------------------------------------------------------------------
            do m=-l,l
              lm=l*(l+1)+m+1
              km=ind_prod(ie,lm,isort)
	        lm0=l+m+1
	        s(ind0+km)=dcmplx(y0(lm0),y1(lm0))
	      enddo
	    enddo
	  enddo   !! over l
	enddo !!  over iatom
	deallocate(y0,y1,t0,s0)
      end
