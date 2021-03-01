      subroutine ws_k_nu_solid(key_u)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      use wannier90_mod      
      implicit none
      include "mpif.h"
      integer, intent(in) :: key_u
      integer :: k,ind_k,npbt1,i_nu,mx2,i,ispin,ind_tau,it,ndd,ind_nu,
     $  iatom_c,iatom,isort,n,j,ii,jj,kk,nn,ierr
      double precision :: ak,dev
      double precision, allocatable :: gx_tau(:)
      complex*16, allocatable :: tmp(:,:),tmp1(:,:)


      mx2=maxel_loc**2
      allocate(w_pr_nu(0:ndiv_r1,0:n_nu))
      w_pr_nu=0.d0
      allocate(w_pr_tau(0:ndiv_r1,0:n_tau/2))
      w_pr_tau=0.d0
      allocate(w_output_int(0:n_nu,3,npnt))
      w_output_int=0.d0
      if(allocated(e_output_int)) deallocate(e_output_int)
      allocate(e_output_int(0:n_nu,4,npnt))
      e_output_int=0.d0
      allocate(e1_output_int(0:n_nu,3,npnt))
      e1_output_int=0.d0
      allocate(p_output_int(0:n_nu,3,npnt))
      p_output_int=0.d0
      allocate(em1_gg0(npnt,0:n_nu))
      em1_gg0=0.d0
      allocate(em1_gg0_low(npnt,0:n_nu))
      em1_gg0_low=0.d0
      npbt1=n_pbmt+nplwgw(1)
      if(allocated(eps_invers_gamma)) deallocate(eps_invers_gamma)
      allocate(eps_invers_gamma(npbt1,npbt1,ndim3_nu))
      eps_invers_gamma=(0.d0,0.d0)
      if(allocated(w_gamma)) deallocate(w_gamma)
      allocate(w_gamma(n_pbmt,n_pbmt,ndim3_nu))
      w_gamma=(0.d0,0.d0)
c     ----------------------------------------------------------------------
C$$$      if(num_orb_proj/=0) then
C$$$        u_loc=0.d0
C$$$      endif
c     -----------------------------------------------------------------------
      allocate(eps_b_mm(4,n_pbtot,ndim3_nu))
      allocate(eps_b_ii(n_pbtot,ndim3_nu))
      eps_b_mm=(0.d0,0.d0)
      eps_b_ii=(0.d0,0.d0)
      dev=0.d0
      e_c_solid_pw=0.d0
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      call timel('**** w at q .ne. 0  start     *******')            
      do ind_k=1,ndim3_k(me3_k+1)
        k=n3_mpi_k(me3_k+1)+ind_k
        ak=dot_product(pnt(:,k),pnt(:,k))
        if(ak>1.d-12) then
          u_rotated_fpb=0.0d0
c$$$          u_at_lpb_rax_k_nu=0.0d0
          w_rotated_fpb=0.0d0
c$$$          w_at_lpb_rax_k_nu=0.0d0
          call ws_k_nu_solid_0(ind_k,key_u)
c$$$          if (myrevision .eq. 1) then
c$$$            call u_w_pade(ind_k, 1)
c$$$            call u_w_pade(ind_k, 2)
c$$$          endif
        endif
      enddo                     !! over k


c     --------- Collect EM1_gg0 ------------------------------------------
      if(nproc_k/=1) then
        call mpi_allreduce_dcmplx
     $    (eps_b_mm,4*n_pbtot*ndim3_nu,mpi_sum,comm_tau)
        call mpi_allreduce_dcmplx
     $    (eps_b_ii,n_pbtot*ndim3_nu,mpi_sum,comm_tau)
        call mpi_allreduce_dble(em1_gg0,npnt*(n_nu+1),mpi_sum,comm_tau)
        call mpi_allreduce_dble
     $    (em1_gg0_low,npnt*(n_nu+1),mpi_sum,comm_tau)
      endif
      if(w_q_zero==0) then
        if(nproc_tau/=1) call mpi_allreduce_dble
     $    (em1_gg0,npnt*(n_nu+1),mpi_sum,comm_k)
        if(nproc_tau/=1) call mpi_allreduce_dble
     $    (em1_gg0_low,npnt*(n_nu+1),mpi_sum,comm_k)
        call output_em1_gg0(em1_gg0,1)
        call output_em1_gg0(em1_gg0_low,0)
      else if(w_q_zero==2) then
        if(nproc_k/=1) then
          call mpi_allreduce_dcmplx
     $      (eps_invers_gamma,npbt1**2*ndim3_nu,mpi_sum,comm_tau)
          call mpi_allreduce_dcmplx
     $      (w_gamma,n_pbmt**2*ndim3_nu,mpi_sum,comm_tau)
        endif
        call eps_invers_gamma_prep2
        if(nproc_tau/=1) call mpi_allreduce_dble
     $    (em1_gg0,npnt*(n_nu+1),mpi_sum,comm_k)
        call output_em1_gg0(em1_gg0)
      endif
c     --------------------------------------------------------------------
      call timel('**** w at q .eq. 0  start    *******')                  
      
      if(me3_k==0)  then
        u_rotated_fpb=0.0d0
c$$$        u_at_lpb_rax_k_nu=0.0d0
        w_rotated_fpb=0.0d0
c$$$        w_at_lpb_rax_k_nu=0.0d0
        call ws_k_nu_solid_0(1,key_u)
c$$$        if (myrevision .eq. 1) then
c$$$          call u_w_pade(1, 1)
c$$$          call u_w_pade(1, 2)
c$$$        endif
      endif


      if (allocated(p_wannier)) deallocate(p_wannier)
c$$$      if((key_u==1) .and.(u_mode .eq. 'band')) deallocate(g_low_tmp)
      deallocate(eps_b_mm,eps_b_ii)
      if(nproc_k/=1) call mpi_allreduce_dble
     $  (e_c_solid_pw,1,mpi_sum,comm_tau)
      if(nproc_tau/=1) call mpi_allreduce_dble
     $  (e_c_solid_pw,1,mpi_sum,comm_k)
c     --- Macroscopic Dielectric matrix -----------------------
      do i_nu=1,n_nu	
        e_output_int(i_nu,2,1)=1.d0/em1_gg0(1,i_nu)
      enddo
      e_output_int(0,2,1)=2*e_output_int(1,2,1)
c     --------------------------------------------------------------     
      if(key_u==1) then

        if(nproc_k/=1) then

        endif
        
      endif
c     -------------------------------------------------------------------
      call output_w
      deallocate(w_pr_nu,w_pr_tau)
      deallocate(em1_gg0)
      deallocate(em1_gg0_low)
      deallocate(w_output_int,p_output_int,e_output_int,e1_output_int)
      if(w_q_zero==2) deallocate(eps_invers_gamma,w_gamma)
      end




