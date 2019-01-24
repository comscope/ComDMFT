	subroutine output_pol_tau_q
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: i,ind_k,k,n_pbt,ind_tau,i_tau,j,i_len
	real*8 :: q
	complex*16 :: wint
	real*8, allocatable :: p_out_int(:,:,:)
	complex*16, allocatable :: work(:,:)
	allocate(p_out_int(0:n_tau/2,3,npnt))
	p_out_int=0.d0
	do ind_k=1,ndim3_k(me3_k+1)
	  k=n3_mpi_k(me3_k+1)+ind_k
	  n_pbt=n_pbmt+nplwgw(k)
	  allocate(work(n_pbt,n_pbt))
	  do ind_tau=1,ndim3_tau
	    i_tau=me3_tau*ndim3_tau+ind_tau-1
	    do j=1,n_pbt
	      do i=1,j-1
	        work(i,j)=dcmplx(p_w_solid_q_tau(i,j,ind_tau,ind_k),
     ^	                       p_w_solid_q_tau(j,i,ind_tau,ind_k))
	        work(j,i)=conjg(work(i,j))
	      enddo
	      work(j,j)=p_w_solid_q_tau(j,j,ind_tau,ind_k)
	    enddo
          call pw_from_pb(ind_k,n_pbt,n_pbt,1,work,wint,3)
	    p_out_int(i_tau,1,k)=real(wint)
	    p_out_int(i_tau,2,k)=real(work(n_pbt,n_pbt))
	    p_out_int(i_tau,3,k)=real(work(1,1))
	  enddo   !! over ind_tau
	  deallocate(work)
	enddo
c --------- Collect W_MT and W_INT  ------------------------------------------
	if(nproc_k/=1) then
	  call dgop(p_out_int,3*npnt*(n_tau/2+1),'  +',comm_pnt)
	endif
	if(nproc_tau/=1) then
	  call dgop(p_out_int,3*npnt*(n_tau/2+1),'  +',comm_pnk)
	endif
	if(maswrk) then
	  i_len=len_trim(allfile)
	  if(ubi==' gw') then
	    open(3,file=allfile(1:i_len)//'_P_tau_q_pw.gw')
	  else if(ubi==' qp') then
	    open(3,file=allfile(1:i_len)//'_P_tau_q_pw.qp')
	  else if(ubi=='gwg') then
	    open(3,file=allfile(1:i_len)//'_P_tau_q_pw.gwg')
	  endif
	  do k=1,npnt
	    q=sqrt(dot_product(pnt(:,k),pnt(:,k)))
	    write(3,*)' Q = ',k,q
	    do i_tau=0,n_tau/2
	      write(3,'(f12.6,3(1x,e15.8))')tau_mesh(i_tau),
     &		     (p_out_int(i_tau,i,k),i=1,3)
	    enddo
	  enddo
	  close(3)
	endif
	deallocate(p_out_int)
      end
