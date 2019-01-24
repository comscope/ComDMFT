      subroutine v_im_r_from_k(ind_tau)
c	ind_tau = 0 - Transform for Coulomb (V) interaction
c	ind_tau > 0 - Transform for Screened Coulomb (W) interaction
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ind_tau
	integer :: k,k0,ind,i0,j0,i,ind_k,ig,j,j1,iatom,isort,ia,ib,
     &	         ind_ir,ir,kst,ind_i
	real*8 :: v(3),fc,pi2
	complex*16 :: cc
	complex*16, allocatable :: v_tmp(:,:,:),v_tmp1(:,:,:,:),tmp1(:,:),
     &	                         tmp2(:,:),tmp3(:,:)
	allocate(v_tmp(nr_full,n_pbmt,maxval(ndim3_k)))
	allocate(tmp1(nr_full,n_pbmt))
	v_tmp=(0.d0,0.d0)
	do ind_k=1,ndim3_k(me3_k+1)
	  k=n3_mpi_k(me3_k+1)+ind_k
	  if(nplwgw(k)/=0) then
	    allocate(tmp2(n_pbmt,nplwgw(k)))
	    if(ind_tau==0) call v_coul_full_mi(nplwgw(k),tmp2,ind_k)
	    if(ind_tau>0) then
	      call w_red_mi(tmp2,p_w_solid_q_tau(1,1,ind_tau,ind_k),
     &	                  nplwgw(k))
          endif
	    tmp1=(0.d0,0.d0)
	    do i=1,nplwgw(k)
	      i0=indpw_gw(i,k)
	      call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
	      tmp1(ind,1:n_pbmt)=tmp2(1:n_pbmt,i)
	    enddo
	    deallocate(tmp2)
	    call fft3(nrdiv(1),nrdiv(2),nrdiv(3),n_pbmt,tmp1,-1)
c ----- Multiplication with the phase factor ------------------
          do i=1,nr_full
            v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &		  	                     +r_pnt(3,i)*r0b0(:,3)
	      fc=-2.d0*pi*dot_product(pnt(:,k),v)
	      cc=exp(dcmplx(0.d0,fc))
	      v_tmp(i,1:n_pbmt,ind_k)=tmp1(i,1:n_pbmt)*cc
	    enddo
	  endif
	enddo   !! over ind_k
	deallocate(tmp1)
c ----------------------------------------------------------------
	allocate(v_tmp1(k_star_max,maxval(ndim3_k),nr_red,n_pbmt))
	pi2=pi+pi
	v_tmp1=(0.d0,0.d0)
	do ind_k=1,ndim3_k(me3_k+1)
	  k0=n3_mpi_k(me3_k+1)+ind_k
     	  do kst=1,k_star(k0)
	    k=k_list(kst,k0)
	    ig=k_group(k)
	    do j=1,nr_red
		    j1=r_pnt_ind(j)
	      j0=ip_r(j1,ig)
		    call sym_vmi_left(ig,v_tmp(j0,:,ind_k),
     &		                  v_tmp1(kst,ind_k,j,:))
	      do iatom=1,natom
	        isort=is(iatom)
	        v=tshift(:,iatom,ig)-tshift_r(:,j1,ig)
	        fc=pi2*dot_product(pnt(:,k),v)
	        cc=dcmplx(cos(fc),sin(fc))
	        ia=iopb(iatom)
	        ib=ia+n_pbmt0(isort)-1
	        v_tmp1(kst,ind_k,j,ia:ib)=v_tmp1(kst,ind_k,j,ia:ib)*cc
	      enddo
	    enddo
	  enddo
	enddo
	deallocate(v_tmp)
c --------------------------------------------------------------------	
	allocate(v_tmp(nqdiv,ndim3_nr_red(me3_k+1),n_pbmt))
	allocate(tmp3(nqdiv,nr_red))
	do j=1,n_pbmt
	  tmp3=(0.d0,0.d0)
	  do ind_k=1,ndim3_k(me3_k+1)
	    k0=n3_mpi_k(me3_k+1)+ind_k
     	    do kst=1,k_star(k0)
	      k=k_list(kst,k0)
	      ig=k_group(k)	    
	      call zone1_number(pnt(:,k),rb0,ndiv,ind)
	      do i=1,nr_red
	        tmp3(ind,i)=v_tmp1(kst,ind_k,i,j)
	      enddo
	    enddo
	  enddo
	  if(nproc_k/=1) call DGOP(tmp3,2*nr_red*nqdiv,'  +',comm_pnt)
	  do ind_i=1,ndim3_nr_red(me3_k+1)
	    i=n3_mpi_nr_red(me3_k+1)+ind_i
	    do k=1,nqdiv
	      v_tmp(k,ind_i,j)=tmp3(k,i)
	    enddo
	  enddo
	enddo
	deallocate(v_tmp1)
	call fft3(ndiv(1),ndiv(2),ndiv(3),n_pbmt*ndim3_nr_red(me3_k+1),
     &	        v_tmp,1)
	v_tmp=v_tmp/dfloat(nqdiv)
c -----------------------------------------------------------------
	if(ind_tau==0) then
	  if(allocated(s_it_r0)) deallocate(s_it_r0)
	  allocate(s_it_r0(nr_red,n_pbmt,ndim3_kk(me3_k+1)))
	endif
	do j=1,n_pbmt
	  tmp3=(0.d0,0.d0)
	  do ind_i=1,ndim3_nr_red(me3_k+1)
	    i=n3_mpi_nr_red(me3_k+1)+ind_i
	    do ir=1,nqdiv
	      tmp3(ir,i)=v_tmp(ir,ind_i,j)
	    enddo
	  enddo
	  if(nproc_k/=1) call DGOP(tmp3,2*nr_red*nqdiv,'  +',comm_pnt)
	  if(ind_tau==0) then
	    do ind_ir=1,ndim3_kk(me3_k+1)
	      ir=n3_mpi_kk(me3_k+1)+ind_ir
	      s_it_r0(:,j,ind_ir)=tmp3(ir,:)
	    enddo
	  else
	    do ind_ir=1,ndim3_kk(me3_k+1)
	      ir=n3_mpi_kk(me3_k+1)+ind_ir
	      w_rs_mi(:,j,ind_ir)=tmp3(ir,:)
	    enddo
	  endif
	enddo	
	deallocate(v_tmp,tmp3)
      end
