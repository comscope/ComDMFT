      subroutine sym_left(ip0,nl,io,u0,s,s0,liml,key)
c	key=1 Y(m)=Sum_m1 D_m_m1 X_m1
c	key=2 Y(m)=Sum_m1 D_m1_m X_m1
	use solid_mod
      implicit none
	integer :: ip0(natom),io(natom),iatom,ind0,jatom,jnd0,nlm,lm,jnd,
     &    	   liml,key,nl,ind
      real*8 :: u0(maxwig)
	complex*16 :: s(*),s0(*)
	real*8, allocatable :: y0(:),t0(:),y1(:)
	allocate(y0(liml),t0(liml),y1(liml))
      do iatom=1,natom     !!  over atoms
	  ind0=io(iatom)-1
	  if(iatom.lt.natom) nlm=io(iatom+1)-io(iatom)
	  if(iatom.eq.natom) nlm=nl-io(iatom)+1
        jatom=ip0(iatom)
	  jnd0=io(jatom)-1
c ------------- Real components ----------------------------------
	  do lm=1,nlm
	    jnd=jnd0+lm
	    t0(lm)=dreal(s0(jnd))
	  enddo
	  call rotate11(t0,y0,nlm,u0,key)
c ------------- Imaginary components ----------------------------------
	  do lm=1,nlm
	    jnd=jnd0+lm
	    t0(lm)=dimag(s0(jnd))
	  enddo
	  call rotate11(t0,y1,nlm,u0,key)
	  do lm=1,nlm
	    ind=ind0+lm
	    s(ind)=dcmplx(y0(lm),y1(lm))
	  enddo
	enddo !!  over iatom
	deallocate(y0,t0,y1)
      end