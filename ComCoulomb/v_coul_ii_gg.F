      subroutine v_coul_ii_gg(ind_k,numpw)

      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ind_k,numpw
      integer :: k,ig,ig0
      double precision :: pi8,pi2a,qg(3),qg2,qg3
      k=n3_mpi_k(me3_k+1)+ind_k
        pi8=8.d0*pi/amega
        pi2a=(pi+pi)/par
        v_ii(:,:,ind_k)=(0.d0,0.d0)
        do ig=1,nplwgw(k)
          ig0=indpw_gw(ig,k)
          qg=(gbs(:,ig0)+pnt(:,k))*pi2a
          qg2=dot_product(qg,qg)
          qg3=qg2
          if(qg2<1.d-6) qg2=1.d0/v_gg_00
          v_ii(ig,ig,ind_k)=pi8/qg2
        enddo                   !! over ig
      end
