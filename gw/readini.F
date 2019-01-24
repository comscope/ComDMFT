      subroutine readini
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use input_mod
      use files_mod
#ifndef MDL_FIXED_FORMAT
      use inp_mod
#endif
#ifdef USE_HDF5
      use hdf5io_mod
#endif
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      character*1  :: chh
      character*20 :: txt
      integer :: iatom,i,j,it,isort,lb,lpb,l,ll,n,ntl(0:10)
      integer :: ng,ngen
      integer :: ierr
      logical :: file_exist


      
      it=0
      allocate(symmop_rotations(3,3,48))
      allocate(symmop_translations(3,48))
      if(maswrk) then
        open(1,file='ini',status='old')
#ifdef MDL_FIXED_FORMAT
        read(1,'(a1)')chh
        read(1,'(19x,i3,10x,i3,10x,i3,9x,i3)')iter_dft,iter_hf,iter_gw,
     &                                        iter_qp
        read(1,'(16x,f5.3,9x,f5.3,8x,f5.3,11x,f5.3)')admix,adspin,
     &                                               adm_gw,acc_it_gw
        read(1,'(16x,i3,11x,f7.4)')iexch,scal_spin
        read(1,'(20x,i4,9x,i4)')nproc_tau,nproc_k
        read(1,'(15x,i1,8x,e11.3,13x,l1,11x,i1)')irel,clight,
     &                                           rel_interst,irel_core
#else
        call inp_hardfail(.true.)
        call inp_errunit(iun)
        call inp_line(1,ierr) ! title line
        call inp_line(1,ierr) ! start control
        call inp_skip(2,ierr)
        call inp_int(iter_dft,ierr)
        call inp_skip(1,ierr)
        call inp_int(iter_hf,ierr)
        call inp_skip(1,ierr)
        call inp_int(iter_gw,ierr)
        call inp_skip(1,ierr)
        call inp_int(iter_qp,ierr)
        call inp_line(1,ierr) ! control, line 2
        call inp_skip(1,ierr)
        call inp_dbl(admix,ierr)
        call inp_skip(1,ierr)
        call inp_dbl(adspin,ierr)
        call inp_skip(1,ierr)
        call inp_dbl(adm_gw,ierr)
        call inp_skip(1,ierr)
        call inp_dbl(acc_it_gw,ierr)
        call inp_line(1,ierr) ! control, line 3
        call inp_skip(1,ierr)
        call inp_int(iexch,ierr)
        call inp_skip(1,ierr)
        call inp_dbl(scal_spin,ierr)
        call inp_line(1,ierr) ! control, line 4
        call inp_skip(1,ierr)
        call inp_int(nproc_tau,ierr)
        call inp_skip(1,ierr)
        call inp_int(nproc_k,ierr)
        call inp_line(1,ierr) ! control, line 5
        call inp_skip(1,ierr)
        call inp_int(irel,ierr)
        call inp_skip(1,ierr)
        call inp_dbl(clight,ierr)
        call inp_skip(1,ierr)
        call inp_log(rel_interst,ierr)
        call inp_skip(1,ierr)
        call inp_int(irel_core,ierr)
#endif
        if(irel==0) irel_core=0
        irel_core=max(irel,irel_core)
        nrel=1
        if(irel.eq.2) nrel=2
        nrel_core=1
        if(irel_core.eq.2) nrel_core=2
#ifdef MDL_FIXED_FORMAT
        read(1,'(22x,f10.2,10x,l)')temperature,restart
        read(1,'(a1)')chh
        read(1,'(10x,a72)')allfile
        read(1,'(11x,a256)')symgen
#else
        call inp_line(1,ierr) ! control, line 6
        call inp_skip(1,ierr)
        call inp_dbl(temperature,ierr)
        call inp_skip(1,ierr)
        call inp_log(restart,ierr)
        call inp_line(1,ierr) ! start files
        call inp_line(1,ierr) ! files, line 2
        call inp_skip(1,ierr)
        call inp_str(allfile,ierr)
        call inp_line(1,ierr," =") ! start sym
        if (inp_len().eq.2) then
c
c         The input line reads: SYM symgen=
c         This happens when the symmetry is P1 and the only generator is
c         the identity, which is implied in all symmetry strings.
c
          symgen=""
        else
          call inp_skip(2,ierr)
          call inp_str(symgen,ierr)
        endif
#endif
        if (symgen.eq.'input') then
#ifdef MDL_FIXED_FORMAT
          read(1,'(34x,i3)')ngroup
#else
          call inp_line(1,ierr) ! start sym
          call inp_skip(4,ierr)
          call inp_int(ngroup,ierr)
#endif
          do ng = 1, ngroup
            read(1,*)
            do i = 1, 3
#ifdef MDL_FIXED_FORMAT
              read(1,'(4x,3(1x,f14.10),6x,f14.10)')
     +             (symmop_rotations(i,j,ng),j=1,3),
     +             symmop_translations(i,ng)
#else
              call inp_line(1,ierr)
              do j = 1, 3
                call inp_dbl(symmop_rotations(i,j,ng),ierr)
              enddo
              call inp_dbl(symmop_translations(i,ng),ierr)
#endif
            enddo
          enddo
        else
c         Symgen holds a string representation of the spacegroup
c         generators. This string will be parsed in subroutine group.
        endif
#ifdef MDL_FIXED_FORMAT
        read(1,'(15x,f11.7,8x,i3,7x,i3,9x,i3)')par,natom,nsort,istruc
#else
        call inp_line(1,ierr) ! start structure
        call inp_skip(2,ierr)
        call inp_dbl(par,ierr)
        call inp_skip(1,ierr)
        call inp_int(natom,ierr)
        call inp_skip(1,ierr)
        call inp_int(nsort,ierr)
        call inp_skip(1,ierr)
        call inp_int(istruc,ierr)
#endif
      endif
      
      if(goparr) then
#ifdef MPI
        call brdcst(msgint,iter_dft,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,iter_hf,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,iter_qp,4,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,admix,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,adspin,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,clight,8,master,MPI_COMM_WORLD)
        call brdcst(msgint,iexch,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,irel,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,irel_core,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,nrel,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,nrel_core,4,master,MPI_COMM_WORLD)
        call brdcstc(allfile,72,master,MPI_COMM_WORLD)
        call brdcst(msgint,nproc_tau,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,nproc_k,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,rel_interst,4,master,MPI_COMM_WORLD)
        call brdcstc(symgen,256,master,MPI_COMM_WORLD)
        call brdcst(msgint,ngroup,4,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,symmop_rotations,8*3*3*48,master,
     &              MPI_COMM_WORLD)
        call brdcst(msgdbl,symmop_translations,8*3*48,master,
     &              MPI_COMM_WORLD)
        call brdcst(msgdbl,par,8,master,MPI_COMM_WORLD)
        call brdcst(msgint,natom,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,nsort,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,istruc,4,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,scal_spin,8,master,MPI_COMM_WORLD)
#endif        
      endif
      allocate(is(natom))
      allocate(tau(3,natom))
      allocate(txtel(nsort))
      allocate(z(nsort))
      allocate(z_dop(nsort))
      allocate(zcor(nsort))
      allocate(magn_shift(nsort))
      allocate(smt(nsort))
      allocate(h(nsort))
      allocate(nrad(nsort))
      allocate(lmb(nsort))
      allocate(lmpb(nsort))
      if(maswrk) then
#ifdef MDL_FIXED_FORMAT
        read(1,'(9x,50i3)')(is(iatom),iatom=1,natom)
        read(1,'(10x,f9.6,6x,f9.6)')b_a,c_a
        read(1,'(8x,3f21.16)')(rbas(i,1),i=1,3)
        read(1,'(8x,3f21.16)')(rbas(i,2),i=1,3)
        read(1,'(8x,3f21.16)')(rbas(i,3),i=1,3)
        do iatom=1,natom
          read(1,'(8x,3f21.16)')(tau(i,iatom),i=1,3)
        enddo
        read(1,'(23x,i4,x,i4,x,i4)')mdiv
        read(1,'(23x,i4,x,i4,x,i4)')nrdiv
        read(1,'(22x,f5.3,14x,f5.3)')cut_lapw_ratio,cut_pb_ratio
        read(1,'(14x,e6.1)')eps_pb
        read(1,'(12x,i4)')nbndf
        read(1,'(14x,f7.3,10x,f7.3,8x,i5)')emindos,emaxdos,ndos
        read(1,'(18x,i4,9x,e6.1)')n_cont_frac,e_small
        read(1,'(10x,l2,25x,l2)')dosplot,bandplot
        read(1,'(14x,3i4,8x,l1,9x,i3,8x,a3)')(ndiv(i),i=1,3),metal,
     &                                           n_k_div,k_line
        read(1,'(14x,f7.3)')v_v0
        read(1,'(17x,f11.6,12x,i5)')b_extval,iter_h_ext
c -------- Transform to internal magnetic field --------------
        read(1,'(14x,3f7.3)')(b_ext(i),i=1,3)
#else
        call inp_line(1,ierr) ! structure, line 2
        call inp_skip(1,ierr)
        do iatom = 1, natom
          call inp_int(is(iatom),ierr)
        enddo
        call inp_line(1,ierr) ! structure, line 3
        call inp_skip(1,ierr)
        call inp_dbl(b_a,ierr)
        call inp_skip(1,ierr)
        call inp_dbl(c_a,ierr)
        do j = 1, 3 ! structure, lines 4-6
          call inp_line(1,ierr) ! structure, line 3
          call inp_skip(1,ierr)
          do i = 1, 3
            call inp_dbl(rbas(i,j),ierr)
          enddo
        enddo
        do iatom = 1, natom ! structure, tau lines
          call inp_line(1,ierr) ! structure, line 3
          call inp_skip(1,ierr)
          do i = 1, 3
            call inp_dbl(tau(i,iatom),ierr)
          enddo
        enddo
        call inp_line(1,ierr) ! start real space meshes
        call inp_skip(4,ierr)
        do i = 1, 3
          call inp_int(mdiv(i),ierr)
        enddo
        call inp_line(1,ierr) ! real space meshes, line 2
        call inp_skip(1,ierr)
        do i = 1, 3
          call inp_int(nrdiv(i),ierr)
        enddo
        call inp_line(1,ierr) ! start basis
        call inp_skip(2,ierr)
        call inp_dbl(cut_lapw_ratio,ierr)
        call inp_skip(1,ierr)
        call inp_dbl(cut_pb_ratio,ierr)
        call inp_line(1,ierr) ! basis, line 2
        call inp_skip(1,ierr)
        call inp_dbl(eps_pb,ierr)
        call inp_line(1,ierr) ! start zones
        call inp_skip(2,ierr)
        call inp_int(nbndf,ierr)
        call inp_line(1,ierr) ! start dos
        call inp_skip(2,ierr)
        call inp_dbl(emindos,ierr)
        call inp_skip(1,ierr)
        call inp_dbl(emaxdos,ierr)
        call inp_skip(1,ierr)
        call inp_int(ndos,ierr)
        call inp_line(1,ierr) ! dos, line 2
        call inp_skip(1,ierr)
        call inp_int(n_cont_frac,ierr)
        call inp_skip(1,ierr)
        call inp_dbl(e_small,ierr)
        call inp_line(1,ierr) ! dos, line 3
        call inp_skip(1,ierr)
        call inp_log(dosplot,ierr)
        call inp_skip(1,ierr)
        call inp_log(bandplot,ierr)
        call inp_line(1,ierr) ! start k_point
        call inp_skip(2,ierr)
        do i = 1, 3
          call inp_int(ndiv(i),ierr)
        enddo
        call inp_skip(1,ierr)
        call inp_log(metal,ierr)
        call inp_skip(1,ierr)
        call inp_int(n_k_div,ierr)
        call inp_skip(1,ierr)
        call inp_str(k_line,ierr)
        call inp_line(1,ierr) ! start multi_scf
        call inp_skip(2,ierr)
        call inp_dbl(v_v0,ierr)
        call inp_line(1,ierr) ! start magnet
        call inp_skip(2,ierr)
        call inp_dbl(b_extval,ierr)
        call inp_skip(1,ierr)
        call inp_int(iter_h_ext,ierr)
        call inp_line(1,ierr) ! magnet, line 2
        call inp_skip(1,ierr)
        do i = 1, 3
          call inp_dbl(b_ext(i),ierr)
        enddo
#endif
        b_ext=b_ext/sqrt(dot_product(b_ext,b_ext))
#ifdef MDL_FIXED_FORMAT
        read(1,'(15x,i5,11x,i5)')n_tau,n_tau_int
        read(1,'(23x,i5,13x,i5,11x,f8.2)')n_omega_exa,n_omega_asy,
     &                                        omega_max
        read(1,'(26x,i2)')interp_omega_d
        read(1,'(17x,i5,10x,i5,8x,f8.2)')n_nu_exa,n_nu_asy,nu_max
        read(1,'(20x,i2)')interp_nu_d
        read(1,'(a1)')chh
#else
        call inp_line(1,ierr) ! start tau mesh
        call inp_skip(3,ierr)
        call inp_int(n_tau,ierr)
        call inp_skip(1,ierr)
        call inp_int(n_tau_int,ierr)
        call inp_line(1,ierr) ! start omega mesh
        call inp_skip(3,ierr)
        call inp_int(n_omega_exa,ierr)
        call inp_skip(1,ierr)
        call inp_int(n_omega_asy,ierr)
        call inp_skip(1,ierr)
        call inp_dbl(omega_max,ierr)
        call inp_line(1,ierr) ! omega mesh, line 2
        call inp_skip(1,ierr)
        call inp_int(interp_omega_d,ierr)
        call inp_line(1,ierr) ! start nu mesh
        call inp_skip(3,ierr)
        call inp_int(n_nu_exa,ierr)
        call inp_skip(1,ierr)
        call inp_int(n_nu_asy,ierr)
        call inp_skip(1,ierr)
        call inp_dbl(nu_max,ierr)
        call inp_line(1,ierr) ! nu mesh, line 2
        call inp_skip(1,ierr)
        call inp_int(interp_nu_d,ierr)
        call inp_line(1,ierr) ! start atomic data
#endif
        maxb=0
        maxpb=0
        maxntle=0
        it=0
        do isort=1,nsort
          read(1,'(a1)')chh
          read(1,'(a1)')chh
#ifdef MDL_FIXED_FORMAT
          read(1,'(6x,i2,7x,i2)')lb,lpb
#else
          call inp_line(1,ierr) ! atomic data, line lmb= lmpb=
          call inp_skip(1,ierr)
          call inp_int(lb,ierr)
          call inp_skip(1,ierr)
          call inp_int(lpb,ierr)
#endif
          read(1,'(a1)')chh
          maxb=max(maxb,lb)
          maxpb=max(maxpb,lpb)
#ifdef MDL_FIXED_FORMAT
          read(1,'(7x,10i3)')(ntl(l),l=0,lb)
#else
          call inp_line(1,ierr) ! atomic data, line ntle=
          call inp_skip(1,ierr)
          do l = 0, lb
            call inp_int(ntl(l),ierr)
          enddo
#endif
          read(1,'(a1)')chh
          it=it+6
          do l=0,lb
            do n=1,ntl(l)
              it=it+1
              read(1,'(a1)')chh
            enddo
            maxntle=max(maxntle,abs(ntl(l)))
          enddo
        enddo
      endif
      if(goparr) then
#ifdef MPI
        call brdcst(msgdbl,cut_lapw_ratio,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,cut_pb_ratio,8,master,MPI_COMM_WORLD)
        call brdcst(msgint,maxb,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,maxpb,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,maxntle,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,metal,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,n_k_div,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,n_cont_frac,4,master,MPI_COMM_WORLD)
        call brdcstc(k_line,3,master,MPI_COMM_WORLD)
        call brdcst(msgint,n_gamma,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,iter_h_ext,4,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,omega_max,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,nu_max,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,e_small,8,master,MPI_COMM_WORLD)
        call brdcst(msgint,dosplot,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,bandplot,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,n_tau,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,n_tau_int,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,interp_omega_d,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,interp_nu_d,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,ndiv,4*3,master,MPI_COMM_WORLD)
        call brdcst(msgint,mdiv,4*3,master,MPI_COMM_WORLD)
        call brdcst(msgint,nrdiv,4*3,master,MPI_COMM_WORLD)
#endif        
      endif
      allocate(konfig(0:maxb,nsort))
      allocate(atoc(maxntle,0:maxb,nsort))
      allocate(idmd(maxntle,0:maxb,nsort))
      allocate(ntle(0:maxb,nsort))
      allocate(augm(maxntle,0:maxb,nsort))
      allocate(correlated(maxntle,0:maxb,nsort))
      allocate(ptnl(maxntle,0:maxb,nsort,2))
      allocate(lim_pb_mt(0:maxpb,nsort))
      h=0.d0
      nrad=0
      z=0.d0
      z_dop=0.d0
      atoc=0.d0
      if(maswrk) then
        do i=1,it
          backspace 1
        enddo
        do isort=1,nsort
#ifdef MDL_FIXED_FORMAT
          read(1,'(8x,a4,4x,f5.1,12x,f6.3)')txtel(isort),z(isort),
     &                                      magn_shift(isort)
          read(1,'(6x,f8.5,3x,f7.4,6x,i5,7x,f6.4)')smt(isort),
     &         h(isort),nrad(isort),z_dop(isort)
          read(1,'(6x,i2,7x,i2)')lmb(isort),lmpb(isort)
          read(1,'(12x,10i3)')(lim_pb_mt(i,isort),i=0,lmpb(isort))
          read(1,'(7x,10i3)')(ntle(l,isort),l=0,lmb(isort))
#else
          call inp_line(1,ierr) ! atomic data, line txtel=
          call inp_skip(1,ierr)
          call inp_str(txtel(isort),ierr)
          call inp_skip(1,ierr)
          call inp_dbl(z(isort),ierr)
          call inp_skip(1,ierr)
          call inp_dbl(magn_shift(isort),ierr)
          call inp_line(1,ierr) ! atomic data, line smt=
          call inp_skip(1,ierr)
          call inp_dbl(smt(isort),ierr)
          call inp_skip(1,ierr)
          call inp_dbl(h(isort),ierr)
          call inp_skip(1,ierr)
          call inp_int(nrad(isort),ierr)
          call inp_skip(1,ierr)
          call inp_dbl(z_dop(isort),ierr)
          call inp_line(1,ierr) ! atomic data, line lmb=
          call inp_skip(1,ierr)
          call inp_int(lmb(isort),ierr)
          call inp_skip(1,ierr)
          call inp_int(lmpb(isort),ierr)
          call inp_line(1,ierr) ! atomic data, line lim_pb_mt=
          call inp_skip(1,ierr)
          do i = 0, lmpb(isort)
            call inp_int(lim_pb_mt(i,isort),ierr)
          enddo
          call inp_line(1,ierr) ! atomic data, line ntle=
          call inp_skip(1,ierr)
          do i = 0, lmb(isort)
            call inp_int(ntle(i,isort),ierr)
          enddo
#endif
          read(1,'(a1)')chh
c ------- Reading the type of the functions for L <= L_b ------
          do l=0,lmb(isort)
            konfig(l,isort)=10
            do n=1,ntle(l,isort)
#ifdef MDL_FIXED_FORMAT
              read(1,'(i3,3x,a3,f7.3,f7.3,5x,a1,4x,i1)')ll,
     &          augm(n,l,isort),atoc(n,l,isort),ptnl(n,l,isort,1),
     &          correlated(n,l,isort),idmd(n,l,isort)
#else
              call inp_line(1,ierr) ! atomic data
              call inp_int(ll,ierr)
              call inp_str(augm(n,l,isort),ierr)
              call inp_dbl(atoc(n,l,isort),ierr)
              call inp_dbl(ptnl(n,l,isort,1),ierr)
              call inp_str(correlated(n,l,isort),ierr)
              call inp_int(idmd(n,l,isort),ierr)
#endif
              ptnl(n,l,isort,2)=ptnl(n,l,isort,1)
              konfig(l,isort)=min(konfig(l,isort),
     &                            int(ptnl(n,l,isort,1)))
            enddo
          enddo
        enddo
        close(1)
      endif
      if(goparr) then
#ifdef MPI
        call brdcst(msgint,n_omega_exa,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,n_omega_asy,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,n_nu_exa,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,n_nu_asy,4,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,omega_max,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,nu_max,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,adm_gw,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,b_a,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,c_a,8,master,MPI_COMM_WORLD)
        call brdcst(msgint,is,4*natom,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,rbas,8*9,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,tau,8*3*natom,master,MPI_COMM_WORLD)
        call brdcst(msgint,nbndf,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,ndos,4,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,emindos,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,emaxdos,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,v_v0,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,b_extval,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,b_ext,8*3,master,MPI_COMM_WORLD)
        call brdcst(msgint,iter_gw,4,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,acc_it_gw,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,temperature,8,master,MPI_COMM_WORLD)
        call brdcst(msgint,restart,4,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,eps_pb,8,master,MPI_COMM_WORLD)
        call brdcstc(txtel,4*nsort,master,MPI_COMM_WORLD)
        call dgop(z,nsort,'  +',MPI_COMM_WORLD)
        call dgop(z_dop,nsort,'  +',MPI_COMM_WORLD)
        call brdcst(msgdbl,magn_shift,8*nsort,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,smt,8*nsort,master,MPI_COMM_WORLD)
        call dgop(h,nsort,'  +',MPI_COMM_WORLD)
        call igop(nrad,nsort,'  +',MPI_COMM_WORLD)
        call dgop(atoc,nsort*(maxb+1)*maxntle,'  +',MPI_COMM_WORLD)
        call brdcstc(correlated,nsort*(maxb+1)*maxntle,master,
     &                  MPI_COMM_WORLD)
        call brdcst(msgint,idmd,4*nsort*(maxb+1)*maxntle,master,
     &                  MPI_COMM_WORLD)
        call brdcst(msgdbl,ptnl,16*nsort*(maxb+1)*maxntle,master,
     &                  MPI_COMM_WORLD)
        call brdcstc(augm,3*nsort*(maxb+1)*maxntle,master,
     &                  MPI_COMM_WORLD)
        call brdcst(msgint,lmb,4*nsort,master,MPI_COMM_WORLD)
        call brdcst(msgint,lmpb,4*nsort,master,MPI_COMM_WORLD)
        call brdcst(msgint,lim_pb_mt,4*nsort*(maxpb+1),master,
     &                  MPI_COMM_WORLD)
        call brdcst(msgint,konfig,4*nsort*(maxb+1),master,
     &                  MPI_COMM_WORLD)
        call brdcst(msgint,ntle,4*nsort*(maxb+1),master,MPI_COMM_WORLD)
#endif        
      endif
      z=z+z_dop
#ifdef USE_HDF5
c      call fopen('gwl.dat',fid_dat)
c      call fopen(trim(allfile)//'.rst',fid_rst)
#endif


      call check_lda_dmft          
c ---------------------------------------------------------------------
      call set_control
      end
