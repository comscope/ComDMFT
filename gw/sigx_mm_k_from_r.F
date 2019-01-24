      subroutine sigx_mm_k_from_r(sig_mm)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	complex*16, intent(in) :: sig_mm(nfun,nfun,ndim3_kk(me3_k+1))
	integer :: ind,ind_k,k,i,j,ind_ir,ir
	complex*16, allocatable :: p_c(:,:)
	allocate(p_c(nqdiv,nfun))
	do j=1,nfun
	  p_c=(0.d0,0.d0)
	  do ind_ir=1,ndim3_kk(me3_k+1)
	    ir=n3_mpi_kk(me3_k+1)+ind_ir
	    do i=1,nfun
	      p_c(ir,i)=sig_mm(i,j,ind_ir)
	    enddo
	  enddo
	  if(nproc_k/=1) call dgop(p_c,2*nqdiv*nfun,'  +',comm_pnt)
	  call fft3(ndiv(1),ndiv(2),ndiv(3),nfun,p_c,-1)
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    call zone1_number(pnt(:,k),rb0,ndiv,ind)
	    do i=1,nfun
	      pvv(i,j,ind_k)=p_c(ind,i)
	    enddo
	  enddo
	enddo
	deallocate(p_c)
      end
