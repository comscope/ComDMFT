      subroutine w_red_mi(wmi,w_v,ind_k,npl)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: ind_k,npl
      double precision, intent(in) :: w_v(nvopt_max,nvopt_max)
      complex*16, intent(out) :: wmi(n_pbmt,npl)
      integer :: k,n,n0
      complex*16, allocatable :: tmp(:,:),tmp1(:,:)
      k=n3_mpi_k(me3_k+1)+ind_k
      n=n_pbmt
      n0=nvopt(k)
      allocate(tmp(n0,n0))
      allocate(tmp1(n,n0))
      call boson_unpack_tau(tmp,n0,w_v,nvopt_max,n0)
      call zgemm('n','n',n,n0,n0,(1.d0,0.d0),qvopt(1,1,ind_k),n_pbtot,
     &  tmp,n0,(0.d0,0.d0),tmp1,n)
      call zgemm('n','c',n,npl,n0,(1.d0,0.d0),tmp1,n,
     &  qvopt(n_pbmt+1,1,ind_k),n_pbtot,(0.d0,0.d0),wmi,n)
      deallocate(tmp,tmp1)
      end
