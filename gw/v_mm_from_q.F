      subroutine v_mm_from_q(v_rs_mm)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      real*8, intent(out) :: v_rs_mm(n_pbmtm,n_pbmtm,ndim_nrr(me3_k+1))
	integer :: k,k0,ind,iatom,isort,i0,j0,ir0,ir,jatom,jsort,ig,
     &	         ind_ir0,ind_i,j,i,ind_j,m0,m,jj0,jj,ind_q,iq,len
      real*8 :: pi2,tt(3),phase
	complex*16, allocatable :: w0(:,:,:),v1(:,:),tr(:),w1(:,:,:)
	pi2=pi+pi
	v_rs_mm=0.d0
	allocate(w0(n_pbmt,n_pbmt,npnt))
	w0=(0.d0,0.d0)
	do ind_q=1,ndim3_k(me3_k+1)
	  iq=n3_mpi_k(me3_k+1)+ind_q
	  call v_coul_full_mt(w0(1,1,iq),ind_q)
      enddo
	if(nproc_k/=1) call dgop(w0,2*npnt*n_pbmt**2,'  +',comm_pnt)
c ------------------------------------------------------------------
      allocate(tr(natom))
	allocate(v1(n_pbmt,n_pbmt))
	allocate(w1(nqdiv,n_pbmt,ndim_pbtt(me3_tau+1)))
      do k=1,nqdiv
        k0=i_kref(k)
	  ig=k_group(k)
	  call zone1_number(pnt(:,k),rb0,ndiv,ind)
	  do iatom=1,natom
	    tt=tshift(:,iatom,ig)
	    phase=pi2*dot_product(pnt(:,k),tt)
	    tr(iatom)=dcmplx(cos(phase),sin(phase))
	  enddo
	  v1=(0.d0,0.d0)
	  do ind_i=1,ndim_pbtt(me3_tau+1)
	    i=n_mpi_pbtt(me3_tau+1)+ind_i
	    v1(:,i)=w0(i,:,k0)
	  enddo
	  call sym_w_mt_side(2,ig,tr,v1(1,n_mpi_pbtt(me3_tau+1)+1),
     &	                   ndim_pbtt(me3_tau+1))
	  if(nproc_tau/=1) call dgop(v1,2*n_pbmt**2,'  +',comm_pnk)
	  v1=transpose(v1)
	  call sym_w_mt_side(1,ig,tr,v1(1,n_mpi_pbtt(me3_tau+1)+1),
     &	                   ndim_pbtt(me3_tau+1))
	  do ind_j=1,ndim_pbtt(me3_tau+1)
	    j=n_mpi_pbtt(me3_tau+1)+ind_j
	    w1(ind,:,ind_j)=v1(:,j)
	  enddo
	enddo   !! over k
	deallocate(tr,v1)
c -----------------------------------------------------------------
	do ind_j=1,ndim_pbtt(me3_tau+1)
	  call fft3(ndiv(1),ndiv(2),ndiv(3),n_pbmt,w1(1,1,ind_j),1)
	enddo
	w1=w1/dfloat(nqdiv)
	m0=n_mpi_pbtt(me3_tau+1)+1
	m=m0+ndim_pbtt(me3_tau+1)-1
	do ind_ir0=1,ndim_nrr(me3_k+1)
	  ir0=n_mpi_nrr(me3_k+1)+ind_ir0
	  ir=nrr_red(1,ir0)
	  iatom=nrr_red(2,ir0)
	  jatom=nrr_red(3,ir0)
	  jsort=is(jatom)
	  j0=iopb(jatom)
	  if(j0>m) cycle
	  j=j0+n_pbmt0(jsort)-1
	  if(j<m0) cycle
	  j0=max(j0,m0)
	  j=min(j,m)
	  jj0=j0-n_mpi_pbtt(me3_tau+1)
	  jj=j-n_mpi_pbtt(me3_tau+1)
	  j0=j0-iopb(jatom)+1
	  j=j-iopb(jatom)+1
	  isort=is(iatom)
	  i0=iopb(iatom)
	  i=i0+n_pbmt0(isort)-1
	  v_rs_mm(1:n_pbmt0(isort),j0:j,ind_ir0)=w1(ir,i0:i,jj0:jj)
	enddo
	len=n_pbmtm**2*ndim_nrr(me3_k+1)
	if(nproc_tau/=1.and.len/=0) call dgop(v_rs_mm,len,'  +',comm_pnk)
	deallocate(w1)
      end
