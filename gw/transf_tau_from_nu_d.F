      subroutine transf_tau_from_nu_d
	use manager_mod
	use parallel_mod
      implicit none
      integer :: i_nu,n0,n1,i,i0,j,ind_tau,i_tau,it,i1_tau,ind_nu
      real*8, allocatable :: tmp(:,:,:),t(:)
      allocate(tmp(0:n_nu_d,0:n_tau/2,2))
      allocate(t(interp_nu_d))
      tmp=0.d0
      do i_nu=0,n_nu
        if(i_nu<=n_nu_exa) then
          n0=i_nu
          n1=i_nu
          t(1)=0.25d0
        else if(i_nu==n_nu) then
          n0=n_nu_d
          n1=n_nu_d
          t(1)=0.25d0
        else
          n0=ind_nu_int(i_nu)
          n1=n0+interp_nu_d-1
          do i=n0,n1
            i0=i-n0+1
            t(i0)=0.25d0
            do j=n0,n1
              if(j==i) cycle
              t(i0)=t(i0)*(w_nu(i_nu)-w_nu_d(j))/(w_nu_d(i)-w_nu_d(j))
            enddo
          enddo
        endif
        do i=n0,n1
          i0=i-n0+1
          do ind_tau=1,ndim3_tau
            i_tau=me3_tau*ndim3_tau+ind_tau-1
            do it=1,2
              tmp(i,i_tau,it)=tmp(i,i_tau,it)
     &                       +t(i0)*tau_from_nu(i_nu,ind_tau,it)
            enddo
          enddo
        enddo
      enddo
      deallocate(t)
	if(nproc_tau/=1) call dgop(tmp,(n_nu_d+1)*(n_tau+2),'  +',
     &	                         comm_pnk)
	do i_tau=0,n_tau/2
	  i1_tau=n_tau-i_tau
	  do ind_nu=1,ndim3_nu_d(me3_tau+1)
	    i_nu=n3_mpi_nu_d(me3_tau+1)+ind_nu-1
	    tau_from_nu_d(ind_nu,i_tau,1)=tmp(i_nu,i_tau,1)
	    tau_from_nu_d(ind_nu,i_tau,2)=tmp(i_nu,i_tau,2)
	    tau_from_nu_d(ind_nu,i1_tau,1)=tmp(i_nu,i_tau,1)
	    tau_from_nu_d(ind_nu,i1_tau,2)=-tmp(i_nu,i_tau,2)
	  enddo
	enddo
      deallocate(tmp)
      end