! If doing dft or hf, set_g00 only calls green_0 and g_mt_from_g_full. 
      subroutine set_g00
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use lda_dmft_mod
      implicit none
      logical :: correlation
      integer :: k,ispin,i,ind_k

      correlation=.false.

      if(ubi==' qp'.or.ubi==' gw') 
     &  correlation=.true.

! If the previous iteration was dft or hf, and the present iteration is qp or gw, then do this.
      if(ubi==' gw'.or.ubi==' qp') then
        if(ubi_0=='dft'.or.ubi_0==' hf') then
          e_qp=e_bnd
          chem_pot_qp=chem_pot
          z_qp=z_bnd
          ev_qp=ev_bnd
          q_qp=(0.d0,0.d0)
          do ispin=1,nspin
            do ind_k=1,ndim3_k(me3_k+1)
              k=n3_mpi_k(me3_k+1)+ind_k
              do i=1,n_bnd(k,ispin)
                q_qp(i,i,ind_k,ispin)=(1.d0,0.d0)
              enddo
            enddo
          enddo
        endif ! (ubi_0=='dft'.or.ubi_0==' hf')
      endif ! (ubi==' gw'.or.ubi==' qp')

! If this is the first gw iteration, do this.
      if(ubi==' gw'.and.ubi_0/=' gw') then
        do ispin=1,nspin
          do ind_k=1,ndim3_k(me3_k+1)
            call set_g_tau_0(ispin,ind_k,0,
     &        g_full(1,1,1,1,ind_k,ispin))
          enddo
        enddo
      endif

! If this is the first qp iteration and the previous iteration was not gw, then do this.
      if(ubi==' qp'.and.ubi_0/=' qp'.and.ubi_0/=' gw') then
        do ispin=1,nspin
          do ind_k=1,ndim3_k(me3_k+1)
            call set_g_tau_0(ispin,ind_k,0,
     &        g_full(1,1,1,1,ind_k,ispin))
          enddo
        enddo
      endif

      if((iter/=1) .or. is_lda_dmft) then

c     ------------- Forming G_FULL(0-) --> G_FULL_0 ------------------------
! green_0 calculates:
! If (dft or hf) and (.not.is_lda_dmft):
!      g_full_00(i,k,ispin)=-green0(de,tau_mesh(n_tau))
! If (dft or hf) and (is_lda_dmft):

!      g_full_0(k) = -green0(k)+deltarhomat1(k)-deltarhomat2(k)
!      g_full_0(k) = -green0(k)+proj_mat(k) * dmft_den_mat(k) * proj_mat^\dagger(k) + (proj_mat(k) * proj_mat^\dagger(k)) * green0 * (proj_mat(k) * proj_mat^\dagger(k))
!      proj_mat(k) = umat(k) * disentangle_mat(k)
!      umat(k) is a rotation basis that is prepared by dft_basis_rotmat_k.  
! If (qp or gw):
!      g_full_0 = an unpacked version of -g_full
! green0 is the Fermi-Dirac function, applied to e_bnd(k)-chem_pot
! g_full_00 and g_full_0 are used by g_mt_from_g_full to create g_loc_0.  It is used by gx_rs_from_ks_mi, gx_rs_from_ks_mm, gx_k_g_r1, and sumi_gw. ro_int_g/sumi_gw uses g_full_00 and g_full_0 to calculate the interstitial sector of the charge density. g_full_00 is also used by seteny.
! -----------------------
! ro_int_g calculates the interstitial sector of the charge density, rointr, and symmetrizes it.
! If doing irel<=1, 
! rointr = \sum_k ev_bnd(k) * g * ev_bnd^\dagger * w(k) / amega, plus symmetrization.
! w(k) = wgt(k)*2.d0/dfloat(nrel*nspin)
! ev_bnd is the interstitial sector of the DFT bands.
! g is g_full_00 if doing dft/hf but not with dmft, and is g_full_0 if doing dmft or qp or gw.
! If doing irel=2, the equations are a bit different.
        call green_0(correlation)
c     ------------- Forming Local Green's Function -------------------------
! g_loc_0(k) = (symm(z_bnd(k) * g(k) * z_bnd^\dagger(k)))^*
! z_bnd is the muffin-tin sector of the dft eigenstates.
! symm means symmetrizing over the crystal group operations.
! Initialize g(k) is g_full_00 if doing dft without dmft, or otherwise g_full_0.
! g_loc_0 is used by get_t, t_t1_x, orb_tot, occ_number, get_pv, spin_orbit_rel. It is also used by the get_t+ro_val pair which calculates the muffin-tin sector of the charge density.
        call g_mt_from_g_full(correlation)

      endif
c     ----------------------------------------------------------------------
      end
