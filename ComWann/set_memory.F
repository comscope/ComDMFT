      subroutine set_memory
      use comwann_mod

      implicit none
      integer :: iatom,isort,
!     $  tmparr(3),
!     $  k0,k1,k2,k3,
!     $  in_out(2,2,2,ndiv(1),ndiv(2),ndiv(3)),
     $  j,j1,j0,nlmb
      character*20 :: ch_atom
!      double precision ::
!     $  tempvec3(3),weig_rc_tmp(ndiv(1),ndiv(2),ndiv(3)),dist_rc_tmp,
!     $  maxdev,gtild(3),v(3)
!      character*300 :: ch,tempchar
!      double precision, allocatable :: tmppnt(:,:),dist_r(:)
      integer, external :: gbs_number

      nlmb=nrel*maxb+1      
      allocate(ffsmt(2,2,maxntle,maxntle,nrel*maxb+1,nsort,nspin))
      allocate(gfun(maxwf,nspin))
      gfun=0.0d0
      allocate(gfund(maxwf,nspin))
      gfund=0.0d0
      allocate(e_bnd(nbndf,npnt,nspin_0))
      
      allocate(n_bnd(npnt,nspin_0))

      allocate(e_qp(nbndf,npnt,nspin_0))
      e_qp=0.0d0
      allocate(ev_bnd(nbasmpw,nbndf,ndim_k(me+1),nspin_0))
      ev_bnd=0.0d0
      allocate(z_bnd(nfun,nbndf,ndim_k(me+1),nspin_0))
      z_bnd=0.0d0
      allocate(ev_qp(nbasmpw,nbndf,ndim_k(me+1),nspin_0))
      ev_qp=0.0d0
      allocate(z_qp(nfun,nbndf,ndim_k(me+1),nspin_0))
      z_qp=0.0d0

      allocate(p_f(0:maxnrad,maxntle,nlmb,nsort,nspin))
      p_f=0.0d0
      allocate(q_f(0:maxnrad,maxntle,nlmb,nsort,nspin))
      q_f=0.0d0
      allocate(v_mt_h(maxmt/nspin_0))
      v_mt_h=0.0d0
      allocate(v_mt_xc(maxmt))
      v_mt_xc=0.0d0      
      allocate(ptnl(maxntle,0:maxb,nsort,2))
      ptnl=0.0d0

      allocate(eny(maxntle,nlmb,nsort,nspin_0))
      eny=0.0d0


!     for wannier input
      allocate(atom_symbols(natom))
      do iatom=1, natom
        isort=is(iatom)
        if (iatom .lt. 10) then
          write(ch_atom,'(i1)') iatom
          atom_symbols(iatom)=trim(txtel(isort))//'_'//
     $      trim(ch_atom)//repeat(" ",16)
          
        elseif (iatom .lt. 100) then
          write(ch_atom,'(i2)') iatom
          atom_symbols(iatom)=trim(txtel(isort))//'_'//
     $      trim(ch_atom)//repeat(" ",15)
        else
          write(ch_atom,'(i3)') iatom
          atom_symbols(iatom)=trim(txtel(isort))//'_'//
     $      trim(ch_atom)//repeat(" ",14)
        endif
      enddo

      allocate(kpt_latt(3,nqdiv))
      
      kpt_latt=0.0d0
!      indtemp=0
      kpt_latt=matmul(transpose(rbas),pnt)

      allocate(t1_x(0:maxb+maxb,
     $  maxlfun,maxlfun,nrel*maxb+1,nsort,nspin))
      t1_x=0.0d0
      if(maxl_cor>=0)
     $  allocate(t_x(nrel*maxb+1,nrel*maxl_cor+1,
     $  0:maxb+maxl_cor,nsort,nspin))
      t_x=0.d00
      
      allocate(g_loc_0(maxel,maxel,natom,nspin_0))
      g_loc_0=0.0d0
      allocate(pcor(maxmtcor))
      pcor=0.0d0
      allocate(qcor(maxmtcor))
      qcor=0.0d0

      allocate(coeff_radial(100,0:maxb,nrel,natom))
      coeff_radial=0.0d0
      allocate(dim_radial(0:maxb,nrel,natom))
      dim_radial=0

      
      allocate(weight_in_froz_win(0:maxb,nrel,nsort))
      weight_in_froz_win=0.0d0
      allocate(weight_low(0:maxb,nrel,nsort))
      weight_low=0.0d0
      allocate(weight_total(0:maxb,nrel,nsort))
      weight_total=0.0d0
      allocate(occ_iwan(0:maxb,nrel,nsort))
      occ_iwan=0.0d0


      allocate(occ_fwan(0:maxb,nrel,nsort))
      occ_fwan=0.0d0

      allocate(weight_fwan(0:maxb,nrel,nsort))
      weight_fwan=0.0d0            

      allocate(bound_radfun(0:maxnrad,2,(nrel*maxb+1),natom,nspin))
      bound_radfun=0.0d0
      allocate(bound_radfun_overlap
     $  (2,maxntle,(nrel*maxb+1),natom,nspin))
      bound_radfun_overlap=0.0d0
      allocate(bound_radfun_coeff
     $  (2,maxntle,(nrel*maxb+1),natom,nspin))
      bound_radfun_coeff=0.0d0

      allocate(sort_degen(nsort))
      sort_degen=0
      
      
      end
      
