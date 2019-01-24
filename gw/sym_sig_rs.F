      subroutine sym_sig_rs(f_red,f)
c     Get F for all R-sites from symmetrized (reduced) representation
c     i.e. to get F(R,tEL,t'E'L') from F_RED(EL,E'L',Rtt') 
c     where Rtt' is reduced to compact number of terms
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	real*8, intent(in) :: f_red(maxel,maxel,nrr_max)
	complex*16, intent(out) :: f(nfun,nfun,ndim3_kk(me3_k+1))
	integer :: j,ig,iatom,isort,jatom,jsort,jnd0,ind0,ir,ir0,l,lp1,ie,
     &	         m,lm,km,lm0,i,in,jn,ind_r
	real*8 :: y0(2*maxb+1),t0(2*maxb+1)
      real*8, allocatable :: tmp(:,:)
      allocate(tmp(maxel,maxel))
      f=(0.d0,0.d0)
	do jatom=1,natom
	  jsort=is(jatom)
	  jnd0=io_lem(jatom)-1
	  do iatom=1,natom
	    isort=is(iatom)
	    ind0=io_lem(iatom)-1
	    do ind_r=1,ndim3_kk(me3_k+1)
	      ir=n3_mpi_kk(me3_k+1)+ind_r
	      ir0=r0_pair(ir,iatom,jatom)
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
	                f(ind0+i,jnd0+km,ind_r)=
     &				      f(ind0+i,jnd0+km,ind_r)+dcmplx(y0(lm0),0.d0)
                    enddo
	            enddo
	          enddo
	        enddo   !! over l
	      enddo  !! over i
	    enddo
	  enddo
	enddo
	deallocate(tmp)
      end