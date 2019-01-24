      subroutine sym_pro_w_side(s0,s,ip0,u0,n)
	use atom_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: ip0(natom),n
	real*8, intent(in) :: s0(n_pbmt,n),u0(maxwig)
	real*8, intent(out) :: s(n_pbmt,n)
	integer :: iatom,ind0,isort,jatom,jnd0,l,ie,lm,km,l2m,iwig
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
	      call dgemm('n','n',l2m,n,l2m,1.d0,u0(iwig),l2m,
     &                 s0(jnd0+km,1),n_pbmt,0.d0,s(ind0+km,1),n_pbmt)
	    enddo
	  enddo   !! over l
      enddo
      end