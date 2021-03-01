      subroutine solid_scf_new

      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      use units_mod

      logical :: we_need_v
      integer :: key_u

      integer :: me3k_ind, ind_k,k,ikstar,kstar,k0,ir1,ir2,ir3,
     $  mr1,mr2,mr3,ir_ind,mr_ind

      we_need_u = .true.
      key_u=1
      we_need_v= .true.
      ubi_num =7
      ubi='dmf'


      if (maswrk) then
        write(iun,*) 'k point information'
        do me3k_ind=1, nproc_k
          do ind_k=1,ndim3_k(me3k_ind)
            k=n3_mpi_k(me3k_ind)+ind_k
            do ikstar=1, k_star(k)
              kstar=k_list(ikstar,k)
              k0=i_kref(kstar)
              write(iun,'(6i5)') me3k_ind, ind_k,ikstar, k,kstar,k0
            enddo
          enddo
        enddo

        do ir1=0,ndiv(1)-1
          do ir2=0,ndiv(2)-1
            do ir3=0,ndiv(3)-1
              mr1=ndiv(1)-ir1
              if (mr1 .eq. ndiv(1)) mr1=0
              mr2=ndiv(2)-ir2
              if (mr2 .eq. ndiv(2)) mr2=0                        
              mr3=ndiv(3)-ir3
              if (mr3 .eq. ndiv(3)) mr3=0                                                
              
              call indexing
     $          (nqdiv,3,ndiv,1,ir_ind,(/ir1+1,ir2+1,ir3+1/))
              call indexing
     $          (nqdiv,3,ndiv,1,mr_ind,(/mr1+1,mr2+1,mr3+1/))
              write(iun,'(a,8i5)') 'fft grid',
     $          ir1,ir2,ir3,mr1,mr2,mr3,ir_ind,mr_ind
            enddo
          enddo
        enddo
      endif

      call prod_bas
      call trans_pb_pw
!              write(100,*)"maxlfun",maxlfun,"maxlfpb",
!     $ maxlfpb,"nsort",nsort,"nspin",nspin,"fi0",fi0
!              call mpi_barrier(mpi_comm_world,mpierr)

      call wan2_pb
!      if (allocated(w_local)) deallocate(w_local)
!      allocate(w_local(nwpmax,nwpmax,ndim3_tau,natom_loc))
!      if (allocated(w_local_nu)) deallocate(w_local_nu)
!      allocate(w_local_nu(nwpmax,nwpmax,ndim3_nu,natom_loc))
c$$$  if (allocated(u_local_nu)) deallocate(u_local_nu)
c$$$  allocate(u_local_nu(nwpmax,nwpmax,ndim3_nu,natom_loc))
!      if (allocated(v_local)) deallocate(v_local)
!      allocate(v_local(nwpmax,nwpmax,natom_loc))
      if(.not.allocated(p_w_solid_q_tau)) then
         allocate(p_w_solid_q_tau(nvopt_max,nvopt_max,ndim3_tau,
     &        ndim3_k(me3_k+1)))
         p_w_solid_q_tau=0.0d0
      endif

!      p_w_solid_q_tau=0.0d0

      if(.not.allocated(p_r_0)) then
         allocate(p_r_0(n_pbtot,n_pbtot,ndim3_tau))
         p_r_0=0.0d0
      endif
      
      call polarizability(key_u)
      call screened_w(key_u)      
      call u_atomic  
      end
