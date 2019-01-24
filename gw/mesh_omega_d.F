      subroutine mesh_omega_d
c --- Here for compatibility when TAU--->OMEGA -------------------------
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
	integer :: i_omega,n
      allocate(ind_omega_d(0:interp_omega_d*(n_omega+1)))	
      allocate(ind_omega_int(n_omega))
      do i_omega=0,n_omega_exa
        ind_omega_d(i_omega)=i_omega
      enddo
      n_omega_d=n_omega_exa
      do i_omega=n_omega_exa+1,n_omega-1
        n=(w_omega(i_omega)*betta_t/pi-1.d0)/2.d0
        if(interp_omega_d==2) then
          if(n==ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n+1
            n_omega_d=n_omega_d+1
          else if(n>ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n
            ind_omega_d(n_omega_d+2)=n+1
            n_omega_d=n_omega_d+2
          endif
        else if(interp_omega_d==3) then
          if(n-1==ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n  ! ref to the Matsubara number
            ind_omega_d(n_omega_d+2)=n+1
            n_omega_d=n_omega_d+2
          else if(n==ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n+1
            n_omega_d=n_omega_d+1
          else if(n-1>ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n-1
            ind_omega_d(n_omega_d+2)=n
            ind_omega_d(n_omega_d+3)=n+1
            n_omega_d=n_omega_d+3
          endif
        else if(interp_omega_d==4) then
          if(n-1==ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n  ! ref to the Matsubara number
            ind_omega_d(n_omega_d+2)=n+1
            ind_omega_d(n_omega_d+3)=n+2
            n_omega_d=n_omega_d+3
          else if(n==ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n+1
            ind_omega_d(n_omega_d+2)=n+2
            n_omega_d=n_omega_d+2
          else if(n+1==ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n+2
            n_omega_d=n_omega_d+1
          else if(n-1>ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n-1
            ind_omega_d(n_omega_d+2)=n
            ind_omega_d(n_omega_d+3)=n+1
            ind_omega_d(n_omega_d+4)=n+2
            n_omega_d=n_omega_d+4
          endif
        else if(interp_omega_d==5) then
          if(n-2==ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n-1  ! ref to the Matsubara number
            ind_omega_d(n_omega_d+2)=n
            ind_omega_d(n_omega_d+3)=n+1
            ind_omega_d(n_omega_d+4)=n+2
            n_omega_d=n_omega_d+4
          else if(n-1==ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n  ! ref to the Matsubara number
            ind_omega_d(n_omega_d+2)=n+1
            ind_omega_d(n_omega_d+3)=n+2
            n_omega_d=n_omega_d+3
          else if(n==ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n+1
            ind_omega_d(n_omega_d+2)=n+2
            n_omega_d=n_omega_d+2
          else if(n+1==ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n+2
            n_omega_d=n_omega_d+1
          else if(n-2>ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n-2
            ind_omega_d(n_omega_d+2)=n-1
            ind_omega_d(n_omega_d+3)=n
            ind_omega_d(n_omega_d+4)=n+1
            ind_omega_d(n_omega_d+5)=n+2
            n_omega_d=n_omega_d+5
          endif
        else if(interp_omega_d==6) then
          if(n-2==ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n-1  ! ref to the Matsubara number
            ind_omega_d(n_omega_d+2)=n
            ind_omega_d(n_omega_d+3)=n+1
            ind_omega_d(n_omega_d+4)=n+2
            ind_omega_d(n_omega_d+5)=n+3
            n_omega_d=n_omega_d+5
          else if(n-1==ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n  ! ref to the Matsubara number
            ind_omega_d(n_omega_d+2)=n+1
            ind_omega_d(n_omega_d+3)=n+2
            ind_omega_d(n_omega_d+4)=n+3
            n_omega_d=n_omega_d+4
          else if(n==ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n+1
            ind_omega_d(n_omega_d+2)=n+2
            ind_omega_d(n_omega_d+3)=n+3
            n_omega_d=n_omega_d+3
          else if(n+1==ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n+2
            ind_omega_d(n_omega_d+2)=n+3
            n_omega_d=n_omega_d+2
          else if(n+2==ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n+3
            n_omega_d=n_omega_d+1
          else if(n-2>ind_omega_d(n_omega_d)) then
            ind_omega_d(n_omega_d+1)=n-2
            ind_omega_d(n_omega_d+2)=n-1
            ind_omega_d(n_omega_d+3)=n
            ind_omega_d(n_omega_d+4)=n+1
            ind_omega_d(n_omega_d+5)=n+2
            ind_omega_d(n_omega_d+6)=n+3
            n_omega_d=n_omega_d+6
          endif
        endif
        ind_omega_int(i_omega)=n_omega_d-interp_omega_d+1 ! gen list
      enddo
      if(ind_omega_d(n_omega_d)<n_omega_maz) then
        n_omega_d=n_omega_d+1
        ind_omega_d(n_omega_d)=n_omega_maz
      endif
      allocate(w_omega_d(0:n_omega_d))
      do i_omega=0,n_omega_d
        n=ind_omega_d(i_omega)
        w_omega_d(i_omega)=(2.d0*n+1.d0)*pi/betta_t
      enddo
	allocate(ndim3_omega_d(nproc_tau))
	allocate(n3_mpi_omega_d(nproc_tau))
	call size_shift_par(n_omega_d+1,nproc_tau,ndim3_omega_d,
     &	                  n3_mpi_omega_d)
	allocate(ndim_omega_d(nproc_k))
	allocate(n_mpi_omega_d(nproc_k))
	call size_shift_par(n_omega_d+1,nproc_k,ndim_omega_d,
     &	                  n_mpi_omega_d)
      call timel('**** MESH_OMEGA_D finished *********')
      end