      subroutine sym_pro_w(s0,s,ip0,u0)
	use atom_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: ip0(natom)
	integer :: iatom,ind0,isort,jatom,jnd0,l,ie,lm,km,l2m,iwig
      real*8 s0(n_pbmt,n_pbmt),s(n_pbmt,n_pbmt),u0(maxwig)
	real*8, allocatable :: pa1(:,:)
	allocate(pa1(n_pbmt,n_pbmt))
	s=0.d0
      do iatom=1,natom     !!  over atoms
	  ind0=iopb(iatom)-1
	  isort=is(iatom)
        jatom=ip0(iatom)
	  jnd0=iopb(jatom)-1
	  do l=0,lmpb(isort)
	    lm=l*l+1   !! first index for the given L
	    l2m=l+l+1
          iwig=l*(2*l-1)*(2*l+1)/3+1
          do ie=1,ntle_pb(l,isort)
            km=ind_prod(ie,lm,isort)
	      call dgemm('n','n',l2m,n_pbmt,l2m,1.d0,u0(iwig),l2m,
     &	               s0(jnd0+km,1),n_pbmt,0.d0,pa1(ind0+km,1),n_pbmt)
	    enddo
	  enddo   !! over l
	enddo
      do iatom=1,natom     !!  over atoms
	  ind0=iopb(iatom)-1
	  isort=is(iatom)
        jatom=ip0(iatom)
	  jnd0=iopb(jatom)-1
	  do l=0,lmpb(isort)
	    lm=l*l+1
	    l2m=l+l+1
          iwig=l*(2*l-1)*(2*l+1)/3+1
          do ie=1,ntle_pb(l,isort)
            km=ind_prod(ie,lm,isort)
	      call dgemm('n','t',n_pbmt,l2m,l2m,1.d0,pa1(1,jnd0+km),
     &		           n_pbmt,u0(iwig),l2m,0.d0,s(1,ind0+km),n_pbmt)
	    enddo
	  enddo   !! over l
	enddo
	deallocate(pa1)
      end
