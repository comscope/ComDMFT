      subroutine sym_pro(ip0,lmx,io,u0,s,s0,nl,liml)
	use solid_mod
      implicit none
	integer :: ip0(natom),lmx(nsort),io(natom),nl,liml,j,
     &	       iatom,ind0,isort,nlm,jatom,jnd0,lm,jnd,i,ind
      real*8 :: u0(maxwig),s(nl,nl),s0(nl,nl)
	real*8, allocatable :: y0(:),t0(:),pa1(:,:)
	allocate(y0(liml),t0(liml),pa1(nl,nl))
	do j=1,nl
        do iatom=1,natom     !!  over atoms
	    ind0=io(iatom)-1
	    isort=is(iatom)
	    nlm=(lmx(isort)+1)**2
          jatom=ip0(iatom)
	    jnd0=io(jatom)-1
	    do lm=1,nlm
	      jnd=jnd0+lm
	      t0(lm)=s0(jnd,j)
	    enddo
	    call rotate11(t0,y0,nlm,u0,1)
	    do lm=1,nlm
	      ind=ind0+lm
	      pa1(ind,j)=y0(lm)
	    enddo
	  enddo !!  over iatom
	enddo
	do i=1,nl
        do iatom=1,natom     !!  over atoms
	    ind0=io(iatom)-1
	    isort=is(iatom)
	    nlm=(lmx(isort)+1)**2
          jatom=ip0(iatom)
	    jnd0=io(jatom)-1
	    do lm=1,nlm
	      jnd=jnd0+lm
	      t0(lm)=pa1(i,jnd)
	    enddo
	    call rotate11(t0,y0,nlm,u0,1)
	    do lm=1,nlm
	      ind=ind0+lm
	      s(i,ind)=s(i,ind)+y0(lm)
	    enddo
	  enddo !!  over ib
	enddo
	deallocate(y0,t0,pa1)
      end