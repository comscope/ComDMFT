      subroutine a_cont_weights(achain,weigd,e,n0)
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n0
      real*8, intent(in) :: e
      complex*16, intent(in) :: achain(0:n_cont_frac,n0,n0,
     &                                 maxval(ndim3_k),
     &                                 nspin)
      real*8, intent(out) :: weigd(n0,npnt,nspin)
	integer :: j,i,ispin,ind_k,k,n,info
	complex*16 :: om0
	integer, allocatable :: ipi(:)
      complex*16, allocatable :: a(:,:),g_om(:,:),g0(:)
      weigd=0.d0
      om0=dcmplx(e,e_small)
      allocate(ipi(n0))
      allocate(g0(n0))
      allocate(a(n0,n0))
      allocate(g_om(n0,n0))
	do ispin=1,nspin
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    n=n_bnd(k,ispin)
c ------ Chain approximation ----------------------------------
		  call pade_fermion_chain_out(n,n0,achain(0,1,1,ind_k,ispin),
     &		                          om0,a,n_cont_frac)
		  g_om=(0.d0,0.d0)
		  do j=1,n
		    do i=1,n
		      a(i,j)=-a(i,j)
		      if(j==i) then
		        a(i,j)=a(i,j)+om0+chem_pot-e_bnd(i,k,ispin)
		        g_om(i,j)=(1.d0,0.d0)
		      endif
		    enddo
		  enddo
          CALL ZGETRF(n,n,a,n0,ipi,INFO)
          CALL ZGETRS('n',n,n,a,n0,ipi,g_om,n0,info)
          do i=1,n
	      weigd(i,k,ispin)=-imag(g_om(i,i))/pi*wgt(k)
	    enddo 
        enddo
      enddo
	deallocate(a,g_om,ipi,g0)
	if(nproc_k/=1) call dgop(weigd,n0*npnt*nspin,'  +',comm_pnt)
	weigd=weigd*2.d0/dfloat(nspin)/nrel
      end
