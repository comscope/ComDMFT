      subroutine etot_x_from_pw(key)
c     key = 1 - SIG_X part
c     key = 2 - P(tau=0) * V part            
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      include "mpif.h"
      integer, intent(in) :: key
      integer :: i,ispin,ind_k,k
      double precision :: ex
      if(key==1) then
        ex=0.d0
        do ispin=1,nspin
          do ind_k=1,ndim3_k(me3_k+1)
            k=n3_mpi_k(me3_k+1)+ind_k
            do i=1,n_bnd(k,ispin)
              ex=ex+wgt(k)*sigx_solid_k(i,i,ind_k,ispin)
            enddo
          enddo
        enddo
        if(nproc_k/=1) call mpi_allreduce_dble(ex,1,mpi_sum,comm_tau)        
        if(nspin==1) ex=2.d0*ex
        e_x_solid_pw=e_x_solid_pw+ex/2.d0
      else if(key==2) then
        ex=0.d0
        if(me3_tau==0) then
          do ind_k=1,ndim3_k(me3_k+1)
            k=n3_mpi_k(me3_k+1)+ind_k
            do i=1,nvopt(k)
              ex=ex+wgt(k)*p_w_solid_q_tau(i,i,1,ind_k)*vopt(i,ind_k)
            enddo
          enddo
        endif
        if(nproc_k/=1) call mpi_allreduce_dble(ex,1,mpi_sum,comm_tau)
        if(nproc_tau/=1) call mpi_allreduce_dble(ex,1,mpi_sum,comm_k)        
        e_x_solid_pw=-ex/4.d0
      endif
      end
