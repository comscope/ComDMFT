      subroutine sym_g_iat_jat_rr(f_red,f,iatom,jatom,ir)
c     Gives G for given iatom_jatom_R 
c     from symmetrized (reduced) representation
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: iatom,jatom,ir
	real*8, intent(in) :: f_red(maxel,maxel,nrr_max)
	real*8, intent(out) :: f(maxel,maxel)
	integer :: j,ig,isort,jsort,ir0,iatom0,jatom0,l,lp1,ie,m,lm,km,
     &	         lm0,i,in,jn
	real*8 :: y0(2*maxb+1),t0(2*maxb+1)
      real*8, allocatable :: tmp(:,:)
      allocate(tmp(maxel,maxel))
      f=0.d0
	isort=is(iatom)
	jsort=is(jatom)
	ir0=r0_pair(ir,iatom,jatom)
	iatom0=nrr_red(2,ir0)
	jatom0=nrr_red(3,ir0)
	ig=ig_pair(ir,iatom,jatom)
	do j=1,lfunm(jsort)
	  do l=0,lmb(isort)
          lp1=l+1
          do ie=1,ntle(l,isort)
	      in=1
	      if(augm(ie,l,isort)/='LOC') in=2
	      do jn=1,in
              do m=-l,l
                lm=l*(l+1)+m+1
                km=indbasa(jn,ie,lm,isort)
	          lm0=l+m+1
	          t0(lm0)=f_red(km,j,ir0)
	        enddo
	        call rotate1(t0,y0,l,u(1,ig),1)
              do m=-l,l
                lm=l*(l+1)+m+1
                km=indbasa(jn,ie,lm,isort)
	          lm0=l+m+1
	          tmp(km,j)=y0(lm0)
	        enddo
	      enddo
	    enddo
	  enddo   !! over l
	enddo   !! over j
	do i=1,lfunm(isort)
	  do l=0,lmb(jsort)
          lp1=l+1
          do ie=1,ntle(l,jsort)
	      in=1
	      if(augm(ie,l,jsort)/='LOC') in=2
	      do jn=1,in
              do m=-l,l
                lm=l*(l+1)+m+1
                km=indbasa(jn,ie,lm,jsort)
	          lm0=l+m+1
	          t0(lm0)=tmp(i,km)
	        enddo
	        call rotate1(t0,y0,l,u(1,ig),1)
              do m=-l,l
                lm=l*(l+1)+m+1
                km=indbasa(jn,ie,lm,jsort)
	          lm0=l+m+1
	          f(i,km)=f(i,km)+y0(lm0)
              enddo
	      enddo
	    enddo
	  enddo   !! over l
	enddo  !! over i
	deallocate(tmp)
      end