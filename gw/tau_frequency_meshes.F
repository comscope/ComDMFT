      subroutine tau_frequency_meshes
c     Create TAU and frequency meshes -------------
      use manager_mod
      use parallel_mod
      use units_mod
      implicit none
      integer :: ind_tau,i_tau,i_omega,i_nu
	real*8 :: h1,h12,h2
	real*8, allocatable :: tran(:,:,:)
c ----------------- TAU mesh --------------------------------------- 
	allocate(tau_mesh(0:n_tau))
	call mesh_tau_cheb
c --------- Here for derivatives at tau=beta ------------------------	
	allocate(df_dtau_0(0:2))
      h1=tau_mesh(1)
      h12=tau_mesh(2)
      h2=h12-h1
      df_dtau_0(0)=-(h1+h12)/h1/h12
      df_dtau_0(1)=h12/h1/h2
      df_dtau_0(2)=-h1/h2/h12
c ----------------- FREQUENCY meshes --------------------------------
      call mesh_nu
c ------------------------------------------------------------------      
      call mesh_omega
c --- Basic AB transformations TAU ---> OMEGA -----------------------
      allocate(ab_omega_from_tau(ndim3_tau,0:n_omega_d,2))
      ab_omega_from_tau=0.d0
      allocate(tran(n_tau/2,1+n_omega_d,2))
      call transf_ferm_from_tau(tran,n_omega_d+1,n_omega_d+1,
     &                          w_omega_d(0))
      do ind_tau=1,ndim3_tau
        i_tau=me3_tau*ndim3_tau+ind_tau-1
        if(i_tau/=0) then
          do i_omega=0,n_omega_d
            if(i_tau/=n_tau/2) ab_omega_from_tau(ind_tau,i_omega,1)=
     &                                        tran(i_tau,i_omega+1,1)
            ab_omega_from_tau(ind_tau,i_omega,2)=
     &                                        tran(i_tau,i_omega+1,2)
          enddo
        endif
      enddo
      deallocate(tran)
c --- Prepare General transformations TAU ---> NU -----------------------
      allocate(ab_nu_from_tau(ndim3_tau,0:n_nu_d,2))
      ab_nu_from_tau=0.d0
      allocate(tran(n_tau/2,1+n_nu_d,2))
      call transf_bos_from_tau(tran,n_nu_d+1,n_nu_d+1,w_nu_d(0))
      do ind_tau=1,ndim3_tau
        i_tau=me3_tau*ndim3_tau+ind_tau-1
        if(i_tau/=0) then
          do i_nu=0,n_nu_d
            ab_nu_from_tau(ind_tau,i_nu,1)=tran(i_tau,i_nu+1,1)
            if(i_tau/=n_tau/2) ab_nu_from_tau(ind_tau,i_nu,2)=
     &                                        tran(i_tau,i_nu+1,2)
          enddo
        endif
      enddo
      deallocate(tran)
      call timel('**** TAU_FREQUNCY_MESHES finished **')
      end