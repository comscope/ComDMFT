      subroutine sym_mp_mt(ig,iq,s,s0,key)
c     Acts like sym_w_left but with reverse t<--->t0 order      
c	key=1 Y(m)=Sum_m1 D_m_m1 X_m1
c	key=2 Y(m)=Sum_m1 D_m1_m X_m1
	use atom_mod
	use solid_mod
	use units_mod
      implicit none
	integer,intent(in) :: ig,key,iq
	complex*16, intent(in) :: s0(n_pbmt)
	complex*16, intent(out) :: s(n_pbmt)
	integer :: iatom,ind0,isort,jatom,jnd0,l,ie,m,lm,km,lm0,i0,i
	real*8 :: pi2,phase
	complex*16 :: cc
	real*8, allocatable :: y0(:),y1(:),t0(:)
	allocate(y0(2*maxpb+1),t0(2*maxpb+1),y1(2*maxpb+1))
	pi2=pi+pi
      do iatom=1,natom     !!  over atoms
	  ind0=iopb(iatom)-1
	  isort=is(iatom)
        jatom=ip(iatom,ig)
	  jnd0=iopb(jatom)-1
	  do l=0,lmpb(isort)
          do ie=1,ntle_pb(l,isort)
c ------------- Real components ----------------------------------
            do m=-l,l
              lm=l*(l+1)+m+1
              km=ind_prod(ie,lm,isort)
	        lm0=l+m+1
	        t0(lm0)=dreal(s0(ind0+km))
	      enddo
	      call rotate1(t0,y0,l,u(1,ig),key)
c ------------- Imaginary components ----------------------------------
            do m=-l,l
              lm=l*(l+1)+m+1
              km=ind_prod(ie,lm,isort)
	        lm0=l+m+1
	        t0(lm0)=dimag(s0(ind0+km))
	      enddo
	      call rotate1(t0,y1,l,u(1,ig),key)
c -------------------------------------------------------------------
            do m=-l,l
              lm=l*(l+1)+m+1
              km=ind_prod(ie,lm,isort)
	        lm0=l+m+1
	        s(jnd0+km)=dcmplx(y0(lm0),y1(lm0))
	      enddo
	    enddo
	  enddo   !! over l
	  i0=iopb(jatom)
	  i=i0+n_pbmt0(isort)-1
	  phase=pi2*dot_product(pnt(:,iq),tshift(:,iatom,ig))
	  cc=dcmplx(cos(phase),sin(phase))
	  s(i0:i)=cc*s(i0:i)
	enddo !!  over iatom
	deallocate(y0,y1,t0)
      end