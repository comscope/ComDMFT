      subroutine p_mm_q_from_r(p_mm,ind_tau,p)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ind_tau
	real*8, intent(in) :: p_mm(n_pbmt,n_pbmt,ndim3_kk(me3_k+1))
      real*8, intent(inout) :: p(n_pbtot,n_pbtot,ndim3_tau,
     &                           maxval(ndim3_k))
	integer :: ind,ind_k,k,i,j,ind_ir,ir
	complex*16, allocatable :: p_c(:,:)
	allocate(p_c(nqdiv,n_pbmt))
	do j=1,n_pbmt
	  p_c=(0.d0,0.d0)
	  do ind_ir=1,ndim3_kk(me3_k+1)
	    ir=n3_mpi_kk(me3_k+1)+ind_ir
	    do i=1,j
	      p_c(ir,i)=p_mm(i,j,ind_ir)
	    enddo
	  enddo
	  if(nproc_k/=1) call dgop(p_c,2*nqdiv*j,'  +',comm_pnt)
	  call fft3(ndiv(1),ndiv(2),ndiv(3),j,p_c,-1)
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    call zone1_number(pnt(:,k),rb0,ndiv,ind)
          do i=1,j-1
	      p(i,j,ind_tau,ind_k)=p(i,j,ind_tau,ind_k)+real(p_c(ind,i))
	      p(j,i,ind_tau,ind_k)=p(j,i,ind_tau,ind_k)+imag(p_c(ind,i))
          enddo
	    p(j,j,ind_tau,ind_k)=p(j,j,ind_tau,ind_k)+real(p_c(ind,j))
        enddo  !! over ind_k
	enddo  !! over j
	deallocate(p_c)
      end
