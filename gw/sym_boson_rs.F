      subroutine sym_boson_rs(f_red,f)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	real*8, intent(in) :: f_red(n_pbmtm,n_pbmtm,nrr_max)
	real*8, intent(out) :: f(n_pbmt,n_pbmt,ndim3_kk(me3_k+1))
	integer :: j,ig,iatom,isort,ind_ir,jatom,jsort,ir,ir0,
     &	         l,lp1,ie,m,lm,km,lm0,i,ip0,jp0,ip00
	real*8 :: y0(2*maxpb+1),t0(2*maxpb+1)
      real*8, allocatable :: tmp(:,:)
      allocate(tmp(n_pbmtm,n_pbmtm))
	f=0.d0
	do jatom=1,natom
	  jsort=is(jatom)
	  jp0=iopb(jatom)-1
	  do iatom=1,natom
	    isort=is(iatom)
	    ip0=iopb(iatom)-1
	    do ind_ir=1,ndim3_kk(me3_k+1)
	      ir=n3_mpi_kk(me3_k+1)+ind_ir
	      ir0=r0_pair(ir,iatom,jatom)
	      ig=ig_pair(ir,iatom,jatom)
	      do j=1,n_pbmt0(jsort)
	        do l=0,lmpb(isort)
                lp1=l+1
                do ie=1,ntle_pb(l,isort)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    km=ind_prod(ie,lm,isort)
	              lm0=l+m+1
	              t0(lm0)=f_red(km,j,ir0)
	            enddo
	            call rotate1(t0,y0,l,u(1,ig),1)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    km=ind_prod(ie,lm,isort)
	              lm0=l+m+1
	              tmp(km,j)=y0(lm0)
	            enddo
	          enddo
	        enddo   !! over l
	      enddo   !! over j
	      do i=1,n_pbmt0(isort)
	        ip00=ip0+i
	        do l=0,lmpb(jsort)
                lp1=l+1
                do ie=1,ntle_pb(l,jsort)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    km=ind_prod(ie,lm,jsort)
	              lm0=l+m+1
	              t0(lm0)=tmp(i,km)
	            enddo
	            call rotate1(t0,y0,l,u(1,ig),1)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    km=jp0+ind_prod(ie,lm,jsort)
	              lm0=l+m+1
	              f(ip00,km,ind_ir)=f(ip00,km,ind_ir)+y0(lm0)
	            enddo
	          enddo
	        enddo   !! over l
	      enddo  !! over i
	    enddo    !! over ir
	  enddo
	enddo
	deallocate(tmp)
      end