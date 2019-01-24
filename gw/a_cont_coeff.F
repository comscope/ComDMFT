      subroutine a_cont_coeff(achain,n0)
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n0
      complex*16, intent(out) :: achain(0:n_cont_frac,n0,n0,
     &                                  maxval(ndim3_k),nspin)
	integer :: i_omega,l,j,i,ii,ispin,ind_k,k,n,ind_omega
	complex*16 :: om0
	integer, allocatable :: nom_chk(:)
	real*8, allocatable :: check_pade(:),s0(:,:,:,:),s0_cheb(:,:,:,:)
      complex*16, allocatable :: sig_om(:,:,:),a(:,:)
      allocate(check_pade(7))
      allocate(nom_chk(7))
      do ii=1,7
        if(ii==1) nom_chk(ii)=0
        if(ii==2) nom_chk(ii)=1
        if(ii==3) nom_chk(ii)=2
        if(ii==4) nom_chk(ii)=n_omega_exa-1
        if(ii==5) nom_chk(ii)=n_omega_exa
        if(ii==6) nom_chk(ii)=n_omega-1
        if(ii==7) nom_chk(ii)=n_omega
      enddo
	allocate(sig_om(n0,n0,0:n_omega))
	allocate(a(n0,n0))
      check_pade=0.d0
	do ispin=1,nspin
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    n=n_bnd(k,ispin)
          sig_om=(0.d0,0.d0)
	    do ind_omega=1,ndim3_omega
	      i_omega=me3_tau*ndim3_omega+ind_omega-1
		    call ferm_unpack_omega(sig_om(1,1,i_omega),
     &	                           sig_c_omega(1,1,1,ind_omega,ind_k,
     &                                         ispin),n,n0,n0)
          enddo
	    if(nproc_tau/=1) call dgop(sig_om,2*n0**2*(n_omega+1),'  +',
     &	                             comm_pnk)
	    if(n_cont_frac<=n_omega_exa) then
c --- Another variant of AC: Continued Fractions approximation ------
            call pade_fermion_chain(sig_om,n,n0,
     &                              achain(0,1,1,ind_k,ispin),
     &                              n_cont_frac)
          else
	      allocate(s0(n,n,2,0:n_omega))
	      allocate(s0_cheb(n,n,2,0:n_omega_asy-2))
	      s0=0.d0
	      do ind_omega=1,ndim3_omega
	        i_omega=me3_tau*ndim3_omega+ind_omega-1
		      s0(:,:,:,i_omega)=
     &		    sig_c_omega(1:n,1:n,:,ind_omega,ind_k,ispin)
            enddo
	      if(nproc_tau/=1) call dgop(s0,2*n*n*(n_omega+1),'  +',
     &	                               comm_pnk)
	      call dgemm('n','t',2*n*n,n_omega_asy-1,n_omega_asy-1,1.d0,
     &                 s0(1,1,1,n_omega_exa+1),2*n*n,cheb_omega(0,1),
     &                 n_omega_asy-1,0.d0,s0_cheb(1,1,1,0),2*n*n)
            call ac_fermion_chain(s0,s0_cheb,n,n0,
     &                            achain(0,1,1,ind_k,ispin),
     &                            n_cont_frac)
            deallocate(s0,s0_cheb)
          endif
c --- Here we check Chain approximation  -----------------------
          do ii=1,7
            l=nom_chk(ii)
            om0=dcmplx(0.d0,w_omega(l))
		    call pade_fermion_chain_out(n,n0,
     &		                            achain(0,1,1,ind_k,ispin),om0,a,
     &                                  n_cont_frac)
		    do j=1,n
		      do i=1,n
		        check_pade(ii)=max(check_pade(ii),
     &		                       abs(sig_om(i,j,l)-a(i,j)))
              enddo
            enddo
          enddo
        enddo   !! over ind_k
      enddo   !! over ispin
	if(nproc_k/=1) call dgop(check_pade,7,'max',comm_pnt)
      if(maswrk) then
        do ii=1,7
          l=nom_chk(ii)
          write(iun,'(a33,i8,a4,e12.5)')
     &       ' Max error in A.Cont for Omega = ',mats_num_omega(l),
     &       ' is ',check_pade(ii)
        enddo
      endif
	deallocate(a,sig_om,nom_chk,check_pade)
      end
