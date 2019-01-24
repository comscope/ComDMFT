      subroutine sym_mm_loc_boson(wmm)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	real*8, intent(inout) :: wmm(n_pbmtm,n_pbmtm,natom)
	integer :: iatom,isort,ig,jatom,j,l,km,lm0,i,ie,lm
	real*8, allocatable :: y0(:),t0(:),tmp(:,:,:),pa1(:,:)
	allocate(tmp(n_pbmtm,n_pbmtm,natom),pa1(n_pbmtm,n_pbmtm),
     & 	       y0(2*maxpb+1),t0(2*maxpb+1))
	tmp=wmm
	wmm=0.d0
	do iatom=1,natom
	  isort=is(iatom)
        do ig=1,ngroup
          jatom=ip(iatom,ig)
	    do j=1,n_pbmt0(isort)
	      do l=0,lmpb(isort)
	        do ie=1,ntle_pb(l,isort)
                do lm0=1,l+l+1
                  lm=l*l+lm0
                  km=ind_prod(ie,lm,isort)
	            t0(lm0)=tmp(km,j,jatom)
	          enddo
	          call rotate1(t0,y0,l,u(1,ig),1)
                do lm0=1,l+l+1
                  lm=l*l+lm0
                  km=ind_prod(ie,lm,isort)
	            pa1(km,j)=y0(lm0)
	          enddo
	        enddo   !! over l
	      enddo
	    enddo    !! over j
	    do i=1,n_pbmt0(isort)
	      do l=0,lmpb(isort)
	        do ie=1,ntle_pb(l,isort)
                do lm0=1,l+l+1
                  lm=l*l+lm0
                  km=ind_prod(ie,lm,isort)
	            t0(lm0)=pa1(i,km)
	          enddo
	          call rotate1(t0,y0,l,u(1,ig),1)
                do lm0=1,l+l+1
                  lm=l*l+lm0
                  km=ind_prod(ie,lm,isort)
	  		  wmm(i,km,iatom)=wmm(i,km,iatom)+y0(lm0)
	          enddo
	        enddo
	      enddo
	    enddo  !! over i
	  enddo  !! over ig
	enddo  !! over iatom_c
	wmm=wmm/ngroup
	deallocate(tmp,pa1,y0,t0)
      end