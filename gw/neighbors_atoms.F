      subroutine neighbors_atoms
c     Creates list of neighbors for each atom (correlated) in Unit cell      
	use atom_mod
	use manager_mod
      use sitegen_mod
	use solid_mod
      implicit none
	integer :: iatom,jatom,i,n,nn,nout,it,lat
	real*8 :: v(3)
      integer, allocatable :: iref(:),ia(:),iv(:)
	real*8, allocatable :: a(:),b(:)
	nn_shells=1
	allocate(nn_sites(natom))
	n=n_rr_max*natom
	allocate(a(n),b(n),iref(n),ia(n),iv(n))
	nn_sites=0
	do it=1,2
	  if(it==2) then
	    max_sites=maxval(nn_sites)
	    allocate(atom_nn(0:max_sites,natom))
	    allocate(lat_nn(0:max_sites,natom))
	  endif
	  do iatom=1,natom
	    i=0
	    do lat=1,n_rr_max
	      do jatom=1,natom
	        i=i+1
	        v=rbs(:,lat)+tau(:,jatom)-tau(:,iatom)
	        a(i)=sqrt(dot_product(v,v))
	        ia(i)=jatom
	        iv(i)=lat
	      enddo
	    enddo
	    call xorder(n,a,nout,b,iref)
	    if(it==1) then
	      nn=0
	      do i=1,n
	        if(iref(i)>1.and.iref(i)<=nn_shells+1) nn=nn+1
	      enddo
	      nn_sites(iatom)=nn
	    else if(it==2) then
	      atom_nn(0,iatom)=iatom
	      lat_nn(0,iatom)=1
	      nn=0
	      do i=1,n
	        if(iref(i)>1.and.iref(i)<=nn_shells+1) then
	          nn=nn+1
	          atom_nn(nn,iatom)=ia(i)
	          lat_nn(nn,iatom)=iv(i)
	        endif
	      enddo
	    endif
	  enddo   !! over iatom
	enddo  !! over it
      deallocate(a,b,iref,ia,iv)
      end
