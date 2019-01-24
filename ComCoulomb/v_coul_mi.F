      subroutine v_coul_mi_old(ind_k)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ind_k
      integer :: k,igg,i,ig
      double precision :: qg(3),pi8,sqpi8,pi2,pi2a,qg2,qq_2,qg3

c     -------------------------------------------------------
      pi8=8.d0*pi
      sqpi8=sqrt(pi8)
      pi2=pi+pi
      pi2a=pi2/par
      v_mi(:,:,ind_k)=(0.d0,0.d0)
      k=n3_mpi_k(me3_k+1)+ind_k
      do ig=1,nplwgw(k)
        igg=indpw_gw(ig,k)
        do i=1,3
          qg(i)=(gbs(i,igg)+pnt(i,k))*pi2a
        enddo
        qg2=dot_product(qg,qg)
        qg3=qg2
        call get_cfgt_old(qg2,qg,ig,v_mi(1,1,ind_k))
        if(qg2<1.d-6) qg2=1.d0/v_gg_00
        qq_2=pi8/qg2/amega
        v_mi(:,ig,ind_k)=v_mi(:,ig,ind_k)*qq_2
      enddo			!! over ig
      end
      
