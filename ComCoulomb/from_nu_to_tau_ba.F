      subroutine from_nu_to_tau_ba(f_nu,n0,f_tau,n1,n)
c     ------- Transform NU ---> TAU  --------------------------
c     Input:
c     F_nu is Hermitian and stored in R*8 square A form
c     Output:
c     F_tau - hermitian, but is stored in square form = A
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      include "mpif.h"
      integer, intent(in) :: n,n0,n1
      double precision, intent(in) :: f_nu(n0,n0,ndim3_nu)
      double precision, intent(out) :: f_tau(n1,n1,ndim3_tau)
      integer :: nn,i,j,ij,ind_tau,i_nu,nblock,nstep,ind_nu,j0,jm
      double precision, allocatable :: at(:,:),an(:,:)
      nblock=1000
      nstep=max(1,nblock/n)
      nn=n*nstep
      allocate(at(nn,ndim3_tau))
      allocate(an(nn,0:n_nu))
      j0=0
      do while(j0<n)
        jm=min(nstep,n-j0)
c     -------------------------- COS transformation ------------------
        an=0.d0
        do ind_nu=1,ndim3_nu
          i_nu=me3_tau*ndim3_nu+ind_nu-1
          ij=0
          do j=j0+1,j0+jm
            do i=1,n
              ij=ij+1
              an(ij,i_nu)=f_nu(i,j,ind_nu)
            enddo
          enddo
        enddo
        if(nproc_tau/=1) call
     $    mpi_allreduce_dble(an,nn*(n_nu+1),mpi_sum,comm_k)        
        call from_nu_to_tau_cos(nn,an,at)
        do ind_tau=1,ndim3_tau
          ij=0
          do j=j0+1,j0+jm
            do i=1,n
              ij=ij+1
              f_tau(i,j,ind_tau)=at(ij,ind_tau)
            enddo
          enddo
        enddo
        j0=j0+jm
      enddo
      deallocate(at,an)
      end
