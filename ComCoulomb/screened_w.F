      subroutine screened_w(key_u)
      use atom_mod
      use manager_mod
      use solid_mod
      use parallel_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: key_u
      integer :: ind_k,k,n_pbt
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),tmp2(:,:)
      allocate(tmp(nvopt_max,nvopt_max))
      allocate(tmp1(n_pbtot,nvopt_max))
      allocate(tmp2(n_pbtot,n_pbtot))

!      write(100,*)"itest",itest,"n3_mpi_k",n3_mpi_k,"n_pbmt",
!     $ n_pbmt,"n_pbtot",n_pbtot,"nbndf_bnd_0",nbndf_bnd_0,
!     $ "nplwgw",nplwgw,"nvopt",nvopt,"nvopt_max",nvopt_max,
!     $ "p_w_solid_q_tau",p_w_solid_q_tau,"qvopt",qvopt,
!     $ "qvopt_tild",qvopt_tild,"we_have_w",we_have_w
!      call mpi_barrier(mpi_comm_world,mpierr)

      
c      if(itest/=0) write(179,*)' Symmetry check for P:'
      do ind_k=1,ndim3_k(me3_k+1)
        k=n3_mpi_k(me3_k+1)+ind_k
        n_pbt=n_pbmt+nplwgw(k)
        call boson_unpack_tau(tmp,nvopt_max,
     &    p_w_solid_q_tau(1,1,1,ind_k),nvopt_max,
     &    nvopt(k))
        call zgemm('n','n',n_pbt,nvopt(k),nvopt(k),(1.d0,0.d0),
     &    qvopt_tild(1,1,ind_k),n_pbtot,tmp,nvopt_max,
     &    (0.d0,0.d0),tmp1,n_pbtot)
        call zgemm('n','c',n_pbt,n_pbt,nvopt(k),(1.d0,0.d0),
     &    tmp1,n_pbtot,qvopt_tild(1,1,ind_k),n_pbtot,
     &    (0.d0,0.d0),tmp2,n_pbtot)
        call check_w_sym(k,tmp2,n_pbtot,2)
      enddo
c      if(itest/=0) write(179,*)' Symmetry check done for P'
      call ws_k_nu_solid(key_u)

c      if(itest/=0) write(179,*)' Symmetry check for W:'
      do ind_k=1,ndim3_k(me3_k+1)
        k=n3_mpi_k(me3_k+1)+ind_k
        n_pbt=n_pbmt+nplwgw(k)
        call boson_unpack_tau(tmp,nvopt_max,
     &    p_w_solid_q_tau(1,1,1,ind_k),nvopt_max,
     &    nvopt(k))
c        if(itest/=0) write(179,*) ind_k, 1
        call zgemm('n','n',n_pbt,nvopt(k),nvopt(k),(1.d0,0.d0),
     &    qvopt(1,1,ind_k),n_pbtot,tmp,nvopt_max,
     &    (0.d0,0.d0),tmp1,n_pbtot)
c        if(itest/=0) write(179,*) ind_k, 2        
        call zgemm('n','c',n_pbt,n_pbt,nvopt(k),(1.d0,0.d0),
     &    tmp1,n_pbtot,qvopt(1,1,ind_k),n_pbtot,
     &    (0.d0,0.d0),tmp2,n_pbtot)
c        if(itest/=0) write(179,*) ind_k, 3
        call check_w_sym(k,tmp2,n_pbtot,2)
c        if(itest/=0) write(179,*) ind_k, 4
      enddo

c      if(itest/=0) write(179,*)' Symmetry check done for W:'
c$$$      if (key_u .eq. 1) then
c$$$        if(irel<2) then
c$$$          if(num_orb_loc/=0) then
c$$$            call w_local_from_full
c$$$c$$$  call w2_from_full ! subroutine necessary for vertex calculation
c$$$          endif
c$$$        endif
c$$$      endif
!      if(nbndf_bnd_0(2)/=0) then
!        call w_coarse_from_fine
!      endif
      deallocate(tmp,tmp1,tmp2)
      we_have_w=.true.
      call timel('**** SCREENED_W finished ***********')
      end
