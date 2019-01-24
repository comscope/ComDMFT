      subroutine from_nu_to_tau_cos(n,a_nu,a_tau)
c     ------- Transform A(NU) ---> A(TAU)  --------------------------
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: a_nu(n,0:n_nu)
      double precision, intent(out) :: a_tau(n,ndim3_tau)
      call dgemm('n','n',n,ndim3_tau,n_nu+1,1.d0,a_nu(1,0),n,
     &  tau_from_nu(0,1,1),n_nu+1,0.d0,a_tau,n)
      end
