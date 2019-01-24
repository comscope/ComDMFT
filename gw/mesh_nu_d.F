      subroutine mesh_nu_d
c --- Here for compatibility when TAU--->NU -------------------------
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
	integer :: i_nu,n
      allocate(ind_nu_d(0:interp_nu_d*(n_nu+1)))	
      allocate(ind_nu_int(n_nu))
      do i_nu=0,n_nu_exa
        ind_nu_d(i_nu)=i_nu
      enddo
      n_nu_d=n_nu_exa
      do i_nu=n_nu_exa+1,n_nu-1
        n=w_nu(i_nu)*betta_t/2.d0/pi
        if(interp_nu_d==2) then
          if(n==ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n+1
            n_nu_d=n_nu_d+1
          else if(n>ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n
            ind_nu_d(n_nu_d+2)=n+1
            n_nu_d=n_nu_d+2
          endif
        else if(interp_nu_d==3) then
          if(n-1==ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n  ! reference to the Matsubara number
            ind_nu_d(n_nu_d+2)=n+1
            n_nu_d=n_nu_d+2
          else if(n==ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n+1
            n_nu_d=n_nu_d+1
          else if(n-1>ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n-1
            ind_nu_d(n_nu_d+2)=n
            ind_nu_d(n_nu_d+3)=n+1
            n_nu_d=n_nu_d+3
          endif
        else if(interp_nu_d==4) then
          if(n-1==ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n  ! reference to the Matsubara number
            ind_nu_d(n_nu_d+2)=n+1
            ind_nu_d(n_nu_d+3)=n+2
            n_nu_d=n_nu_d+3
          else if(n==ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n+1
            ind_nu_d(n_nu_d+2)=n+2
            n_nu_d=n_nu_d+2
          else if(n+1==ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n+2
            n_nu_d=n_nu_d+1
          else if(n-1>ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n-1
            ind_nu_d(n_nu_d+2)=n
            ind_nu_d(n_nu_d+3)=n+1
            ind_nu_d(n_nu_d+4)=n+2
            n_nu_d=n_nu_d+4
          endif
        else if(interp_nu_d==5) then
          if(n-2==ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n-1  ! reference to the Matsubara number
            ind_nu_d(n_nu_d+2)=n
            ind_nu_d(n_nu_d+3)=n+1
            ind_nu_d(n_nu_d+4)=n+2
            n_nu_d=n_nu_d+4
          else if(n-1==ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n  ! reference to the Matsubara number
            ind_nu_d(n_nu_d+2)=n+1
            ind_nu_d(n_nu_d+3)=n+2
            n_nu_d=n_nu_d+3
          else if(n==ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n+1
            ind_nu_d(n_nu_d+2)=n+2
            n_nu_d=n_nu_d+2
          else if(n+1==ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n+2
            n_nu_d=n_nu_d+1
          else if(n-2>ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n-2
            ind_nu_d(n_nu_d+2)=n-1
            ind_nu_d(n_nu_d+3)=n
            ind_nu_d(n_nu_d+4)=n+1
            ind_nu_d(n_nu_d+5)=n+2
            n_nu_d=n_nu_d+5
          endif
        else if(interp_nu_d==6) then
          if(n-2==ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n-1  ! reference to the Matsubara number
            ind_nu_d(n_nu_d+2)=n
            ind_nu_d(n_nu_d+3)=n+1
            ind_nu_d(n_nu_d+4)=n+2
            ind_nu_d(n_nu_d+5)=n+3
            n_nu_d=n_nu_d+5
          else if(n-1==ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n  ! reference to the Matsubara number
            ind_nu_d(n_nu_d+2)=n+1
            ind_nu_d(n_nu_d+3)=n+2
            ind_nu_d(n_nu_d+4)=n+3
            n_nu_d=n_nu_d+4
          else if(n==ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n+1
            ind_nu_d(n_nu_d+2)=n+2
            ind_nu_d(n_nu_d+3)=n+3
            n_nu_d=n_nu_d+3
          else if(n+1==ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n+2
            ind_nu_d(n_nu_d+2)=n+3
            n_nu_d=n_nu_d+2
          else if(n+2==ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n+3
            n_nu_d=n_nu_d+1
          else if(n-2>ind_nu_d(n_nu_d)) then
            ind_nu_d(n_nu_d+1)=n-2
            ind_nu_d(n_nu_d+2)=n-1
            ind_nu_d(n_nu_d+3)=n
            ind_nu_d(n_nu_d+4)=n+1
            ind_nu_d(n_nu_d+5)=n+2
            ind_nu_d(n_nu_d+6)=n+3
            n_nu_d=n_nu_d+6
          endif
        endif
        ind_nu_int(i_nu)=n_nu_d-interp_nu_d+1 ! ref to the general list
      enddo  !! over i_nu
      if(ind_nu_d(n_nu_d)<n_nu_maz) then
        n_nu_d=n_nu_d+1
        ind_nu_d(n_nu_d)=n_nu_maz
      endif
      allocate(w_nu_d(0:n_nu_d))
      do i_nu=0,n_nu_d
        n=ind_nu_d(i_nu)
        w_nu_d(i_nu)=2.d0*n*pi/betta_t
      enddo
	allocate(ndim3_nu_d(nproc_tau))
	allocate(n3_mpi_nu_d(nproc_tau))
	call size_shift_par(n_nu_d+1,nproc_tau,ndim3_nu_d,n3_mpi_nu_d)
	allocate(ndim_nu_d(nproc_k))
	allocate(n_mpi_nu_d(nproc_k))
	call size_shift_par(n_nu_d+1,nproc_k,ndim_nu_d,n_mpi_nu_d)
      call timel('**** MESH_NU_D finished ************')
      end