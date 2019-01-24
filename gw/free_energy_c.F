      subroutine free_energy_c(ispin,ind_k)
c - Evaluates  - Tr ( Ln G  - Ln Gx )   ------------------
      use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: ispin,ind_k
	integer :: i_omega,i,ind_omega,n,k
	real*8 :: w,et1
	complex*16 :: gx
	real*8, allocatable :: d(:)
	complex*16, allocatable :: p1(:,:),p2(:,:),sigc(:,:)
	k=n3_mpi_k(me3_k+1)+ind_k
	n=n_bnd(k,ispin)
	allocate(p1(n,n))
	allocate(p2(n,n))
	allocate(sigc(n,n))
	allocate(d(n))
	et1=0.d0
	do ind_omega=1,ndim3_omega
	  i_omega=me3_tau*ndim3_omega+ind_omega-1
	  w=w_omega(i_omega)
	  call ferm_unpack_omega(sigc,
     &		                   sig_c_omega(1,1,1,ind_omega,ind_k,ispin),
     &                         n,nbndf,n)
	  do i=1,n
	    gx=(1.d0,0.d0)/(dcmplx(0.d0,w)+chem_pot-e_bnd(i,k,ispin))
	    p1(i,:)=-gx*sigc(i,:)
	  enddo
	  do i=1,n
	    p1(i,i)=p1(i,i)+(1.d0,0.d0)
	  enddo
	  p2=transpose(conjg(p1))
	  call zgemm('n','n',n,n,n,(1.d0,0.d0),p1,n,p2,n,(0.d0,0.d0),sigc,
     & 	           n)
	  if(me_kk_pbr==0) call eig_val_solver(n,n,sigc,d)
	  d=log(d)
        et1=et1+et_sum_sc_gx(i_omega)*sum(d)/2.d0
      enddo   !!! over ind_omega
	if(nproc_tau/=1) call DGOP(et1,1,'  +',comm_pnk)
      d_free_c=d_free_c-wgt(k)*et1*2.d0/nspin/nrel*2.d0
	deallocate(p1,p2,sigc,d)
      end
      
      
      
      
      
      subroutine free_energy_c_0(ispin,ind_k)
c - Evaluates  - Tr ( Ln G  - Ln Gx )   ------------------
      use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: ispin,ind_k
	integer :: i_omega,i,ind_omega,n,nn,k
	real*8 :: w,et,et1,s
	complex*16 :: gx
	complex*16, allocatable :: p1(:,:),p2(:,:),sigc(:,:)
	allocate(p1(nbndf,nbndf))
	allocate(p2(nbndf,nbndf))
	allocate(sigc(nbndf,nbndf))
	k=n3_mpi_k(me3_k+1)+ind_k
	n=n_bnd(k,ispin)
	et1=0.d0
	do ind_omega=1,ndim3_omega
	  i_omega=me3_tau*ndim3_omega+ind_omega-1
	  w=w_omega(i_omega)
	  call ferm_unpack_omega(sigc,
     &		                   sig_c_omega(1,1,1,ind_omega,ind_k,ispin),
     &                         n,nbndf,nbndf)
	  et=0.d0
	  do i=1,n
	    gx=(1.d0,0.d0)/(dcmplx(0.d0,w)+chem_pot-e_bnd(i,k,ispin))
	    p1(i,:)=gx*sigc(i,:)
	  enddo
	  p2=p1
	  do nn=1,1000
	    s=0.d0
	    do i=1,n
	      s=s+p2(i,i)
	    enddo
	    s=s/nn
	    et=et+s
	    if(abs(s)<1.d-8) exit
	    call zgemm('n','n',n,n,n,(1.d0,0.d0),p2,nbndf,p1,nbndf,
     &               (0.d0,0.d0),sigc,nbndf)
	    p2=sigc
	  enddo
        et1=et1+et_sum_sc_gx(i_omega)*et
      enddo   !!! over ind_omega
	if(nproc_tau/=1) call DGOP(et1,1,'  +',comm_pnk)
      d_free_c=d_free_c+wgt(k)*et1*2.d0/nspin/nrel*2.d0
	deallocate(p1,p2,sigc)
      end
