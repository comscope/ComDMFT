      subroutine from_tau_to_nu_cos_d(f_nu,f_tau,n)
c     ------- Transform TAU ---> NU representation ----------------------
c     This is a basic operation A(v) = Int_0^beta/2 A(tau) cos(v*tau)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      include "mpif.h"
      integer, intent(in) :: n
      double precision, intent(in) :: f_tau(n,ndim3_tau)
      double precision, intent(out) :: f_nu(n,0:n_nu_d)
      call dgemm('n','n',n,n_nu_d+1,ndim3_tau,1.d0,f_tau,n,
     &  ab_nu_from_tau(1,0,1),ndim3_tau,0.d0,f_nu,n)
      if(nproc_tau/=1) call mpi_allreduce_dble
     $  (f_nu,n*(n_nu_d+1),mpi_sum,comm_k)      
      end
