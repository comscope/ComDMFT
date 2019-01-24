      subroutine sym_a_2(ev_new,k,ev_bnd0,neig,ig,k0)
      
c     k0 here is always from NPNT set (not npnt)      
      use comwann_mod
      implicit none
      integer, intent(in) :: k,neig,ig,k0
      complex*16, intent(in) :: ev_bnd0(nbasmpw,neig)
      complex*16, intent(out) :: ev_new(nbasmpw,neig)
      integer :: j,ibnd,j0,j1,gbs_number,i1,i2,ii,i0
      double precision :: v(3),pi2,phase,gtild(3)
      complex*16 :: cc, ev_bnd1(nbasmpw,neig)
      if(k==k0) then
        ev_new=ev_bnd0
        return
      else
        ev_new=0.0d0
      endif
      pi2=pi+pi
      call rotate(pnt(1,k0),pnt(2,k0),pnt(3,k0),v(1),v(2),v(3),
     &  u(2,ig),2)
      gtild=pnt(:,k)-v
c     -------- Symmetrization into another q - point ----------------------

      if (irel .le. 1) then
        do j=1,nbask(k0)
          j0=indgb(j,k)         !! G in GBS-list
          v=gbs(:,j0)+gtild     !! G+G_A
          j1=gbs_number(v)      !! G+G_A in GBS-list
          j1=gbs_sym(j1,ig)     !! A^(-1)*(G+G_A) in GBS-list
          j1=iplf_bk(j1,k0)     !! A^(-1)*(G+G_A) in indgb list
          v=pnt(:,k)+gbs(:,j0)  !! k+G
          phase=-pi2*dot_product(v,shift(:,ig))
          cc=dcmplx(cos(phase),sin(phase))
          do ibnd=1,neig
            ev_new(j,ibnd)=cc*ev_bnd0(j1,ibnd)
          enddo
        enddo                   !! over j
      else
        ev_bnd1=0.0d0
        do ii=1, nrel
          i0=(nbask(k0)/nrel)*(ii-1)
          do j=1,nbask(k0)/nrel
            j0=indgb(j,k)       !! G in GBS-list
            v=gbs(:,j0)+gtild   !! G+G_A
            j1=gbs_number(v)    !! G+G_A in GBS-list
            j1=gbs_sym(j1,ig)   !! A^(-1)*(G+G_A) in GBS-list
            j1=iplf_bk(j1,k0)
            v=pnt(:,k)+gbs(:,j0) !! k+G
            phase=-pi2*dot_product(v,shift(:,ig))
            cc=dcmplx(cos(phase),sin(phase))
            do ibnd=1,neig
              ev_bnd1(j+i0,ibnd)=cc*ev_bnd0(j1+i0,ibnd)
            enddo
          enddo                 !! over j        
        enddo
        do j=1,nbask(k0)/nrel        
          i1=j
          i2=nbask(k0)/nrel+j
          do ibnd=1,neig
            ev_new(i1,ibnd)=conjg(uj(1,ig))*ev_bnd1(i1,ibnd)
     &        +conjg(uj(3,ig))*ev_bnd1(i2,ibnd)
            ev_new(i2,ibnd)=conjg(uj(2,ig))*ev_bnd1(i1,ibnd)
     &        +conjg(uj(4,ig))*ev_bnd1(i2,ibnd)
          enddo
        enddo                
      endif
      end

      
