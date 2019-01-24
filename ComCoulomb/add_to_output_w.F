      subroutine add_to_output_w(ind_k,k,n_pbt,w_c_nu,key)
c     key = 0 Polarizability
c     key = 1 W_full - V
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ind_k,k,n_pbt,key
      double precision, intent(in) :: w_c_nu(n_pbt,n_pbt,ndim3_nu)
      integer :: ind_nu,i_nu
      complex*16 :: wint
      complex*16, allocatable :: work(:,:)
      if(k<=npnt) then
        allocate(work(n_pbt,n_pbt))
        do ind_nu=1,ndim3_nu
          i_nu=me3_tau*ndim3_nu+ind_nu-1
          call unpack_hermit(work,w_c_nu(1,1,ind_nu),n_pbt,n_pbt,n_pbt,
     $      (0.d0,0.d0),(1.d0,0.d0))
          if(key==0) then
            call pw_from_pb(ind_k,n_pbt,n_pbt,work,wint,3)
            p_output_int(i_nu,1,k)=dble(wint)
            p_output_int(i_nu,2,k)=dble(work(n_pbt,n_pbt))
            p_output_int(i_nu,3,k)=dble(work(1,1))
          else if(key==1) then
            call pw_from_pb(ind_k,n_pbt,n_pbt,work,wint,1)
            w_output_int(i_nu,1,k)=dble(wint)
            w_output_int(i_nu,2,k)=dble(work(n_pbt,n_pbt))
            w_output_int(i_nu,3,k)=dble(work(1,1))
          endif
        enddo
        deallocate(work)
      endif
      if(key==1) then
c     ---- Contribution for plotting W along a line in real space -----	
c     call w_prn(w_c_nu,pw,n_pbt,k)
      endif
      end
