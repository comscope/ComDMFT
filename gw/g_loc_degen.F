	subroutine g_loc_degen(g,n,ngen,ngen_size,ind_g,ind_g_back)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	complex*16, intent(in) :: g(maxel_loc,maxel_loc,ndim3_omega,nspin)
	integer, intent(out) :: ngen,ngen_size(n),ind_g(n),ind_g_back(n,n)
	integer :: i,j,ispin,ind_omega,i_omega
	real*8 :: s,eps
	integer, allocatable :: iflag(:)
	complex*16, allocatable :: go(:,:,:,:)
	eps=1.d-8
	allocate(go(n,n,0:n_omega,nspin))
	go=(0.d0,0.d0)
	do ispin=1,nspin
	  do ind_omega=1,ndim3_omega
	    i_omega=me3_tau*ndim3_omega+ind_omega-1
	    go(:,:,i_omega,ispin)=g(1:n,1:n,ind_omega,ispin)
	  enddo
	enddo
	if(nproc_tau/=1) call dgop(go,2*n*n*(n_omega+1)*nspin,'  +',
     &	                         comm_pnk)
	allocate(iflag(n))
	iflag=0
	ngen=0
	ngen_size=0
	do j=1,n
	  if(iflag(j)==1) cycle
	  ngen=ngen+1
	  do i=1,n
	    if(iflag(i)==1) cycle
	    s=0.d0
	    do ispin=1,nspin
	      do i_omega=0,n_omega
	        s=max(s,abs(go(j,j,i_omega,ispin)-go(i,i,i_omega,ispin)))
	      enddo
	    enddo
	    if(s<eps) then
	      ind_g(i)=ngen
	      iflag(i)=1
	      ngen_size(ngen)=ngen_size(ngen)+1
	      ind_g_back(ngen_size(ngen),ngen)=i
	    endif
	  enddo
	enddo
	deallocate(iflag,go)
      end