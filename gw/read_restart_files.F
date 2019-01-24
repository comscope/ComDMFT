c-----------------------------------------------------------------------
c
      subroutine read_restart_ubi_rst
c
c     Read the type of the current calculation phase and the current
c     iteration on the restart file. 
c
c     Note that the program goes through the following phases in order:
c     1. DFT           1                   -  iter_dft
c     2. Hartree-Fock  iter_dft+1          -  iter_dft+iter_hf
c     3. GW            iter_dft+iter_hf+1  -  iter_dft+iter_hf+
c                                             iter_gw
c     4. QP            iter_dft+iter_hf+   -  iter_dft+iter_hf+
c                      iter_gw+1              iter_gw+iter_qp
c
      use manager_mod
      use parallel_mod
      use files_mod
#ifdef USE_HDF5
#ifndef MPI
The code expects to be using MPI if you want to use HDF5
#endif
#endif
#ifdef USE_HDF5
      use hdf5io_mod
#endif
      implicit none
#ifdef MPI
      include "mpif.h"
#endif
#ifdef USE_HDF5
      call fopen(trim(allfile)//'.rst',fid_rst)

      call fread(MPI_COMM_WORLD,ubi,"ubi",fid_rst)
      call fread(MPI_COMM_WORLD,iter,"iteration",fid_rst)
      iter = iter + 1
      call fclose(fid_rst)
#else
      if (maswrk) then
        open(3,file=trim(allfile)//'.rst',status='old',
     &       form="unformatted")
        read(3)ubi,iter
        close(3)
        write(iun,*)'read restart: ',ubi,' step ',iter
c
c       The stored iteration number is the last iteration completed.
c       Therefore the restarted calculation starts at the next
c       iteration. However... see also subroutine read_wfn_effhf_rst
c       for the handling of the variables v_h_bb, v_h_bb_old,
c       v_xc_bb, and v_xc_bb_old.
c
        iter = iter + 1
      endif
#ifdef MPI
      if (goparr) then
        call brdcst(msgchr,ubi,3,master,MPI_COMM_WORLD)
        call brdcst(msgint,iter,4,master,MPI_COMM_WORLD)
      endif
#endif
#endif
      ubi_0 = ubi
      end
c
c-----------------------------------------------------------------------
c
      subroutine read_chemical_potential_rst
c
c     Read the chemical potential from the restart files if required
c     for the current calculation phase.
c
      use manager_mod
      use parallel_mod
      use files_mod
#ifdef USE_HDF5
      use hdf5io_mod
#endif
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      character*3 ubi_l
      ubi_l=adjustl(ubi)
c
#ifdef USE_HDF5
      call fopen(trim(allfile)//'.rst',fid_rst)
      call fread(MPI_COMM_WORLD,chem_pot,"chemical_potential/"//
     &           trim(ubi_l)//"/chem_pot",fid_rst)
      if (ubi.eq.' qp') then
        call fread(MPI_COMM_WORLD,chem_pot_qp,"chemical_potential/"//
     &             trim(ubi_l)//"/chem_pot_qp",fid_rst)
      endif
      call fclose(fid_rst)
#else
      if (maswrk) then
        open(3,file=trim(allfile)//'_chemical_potential_'//
     &       trim(ubi_l)//'.rst',status='old',form="unformatted")
        read(3)chem_pot
        if (ubi.eq.' qp') then
          read(3)chem_pot_qp
        endif
        close(3)
      endif
#ifdef MPI
      if (goparr) then
        call brdcst(msgdbl,chem_pot,8,master,MPI_COMM_WORLD)
        if (ubi.eq.' qp') then
          call brdcst(msgdbl,chem_pot_qp,8,master,MPI_COMM_WORLD)
        endif
      endif
#endif

#endif
      end
c
c-----------------------------------------------------------------------
c
      subroutine read_ro_rst
c
c     Read the density of the valence electrons from the restart file
c
      use manager_mod
      use parallel_mod
      use solid_mod
      use atom_mod
      use files_mod
#ifdef USE_HDF5
      use hdf5io_mod
#endif
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      character*3 ubi_l
      integer :: ng(2)
      integer :: nv(2)
      ubi_l=adjustl(ubi)
c
#ifdef USE_HDF5
      call fopen(trim(allfile)//'.rst',fid_rst)
      nv = 0
      ng(1) = size(ro,1)
      call fread(MPI_COMM_WORLD,ro,ng(1),nv(1),"el_density/"//
     &           trim(ubi_l)//"/ro",fid_rst)
      ng(1) = size(rointr,1)
      ng(2) = size(rointr,2)
      call fread(MPI_COMM_WORLD,rointr,ng,nv,"el_density/"//
     &           trim(ubi_l)//"/rointr",fid_rst)
      if (magn.eq.2) then
        ng(1) = size(spmt,1)
        call fread(MPI_COMM_WORLD,spmt,ng(1),nv(1),"el_density/"//
     &             trim(ubi_l)//"/spmt",fid_rst)
        ng(1) = size(spintr,1)
        ng(2) = size(spintr,2)
        call fread(MPI_COMM_WORLD,spintr,ng,nv,"el_density/"//
     &             trim(ubi_l)//"/spintr",fid_rst)
      endif
      call fclose(fid_rst)
#else
      if (maswrk) then
        open(3,file=trim(allfile)//'_el_density_'//trim(ubi_l)//
     &       '.rst',status="old",form="unformatted")
        read(3)ro
        read(3)rointr
        if (magn.eq.2) then
          read(3)spmt
          read(3)spintr
        endif
        close(3)
      endif
c
c     ro     - double precision
c     rointr - double complex
c     spmt   - double precision
c     spintr - double complex
c
#ifdef MPI
      if (goparr) then
        call brdcst(msgdbl,ro,8*maxmt,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,rointr,16*nplwro*nspin_0,master,
     &              MPI_COMM_WORLD)
        if (magn.eq.2) then
          call brdcst(msgdbl,spmt,8*maxmtb,master,MPI_COMM_WORLD)
          call brdcst(msgdbl,spintr,3*16*nplwro,master,
     &                MPI_COMM_WORLD)
        endif
      endif
#endif
#endif
      end
c
c-----------------------------------------------------------------------
c
      subroutine read_ro_core_rst
c
c     Read the core electron densities from the restart files
c
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use files_mod
#ifdef USE_HDF5
      use hdf5io_mod
#endif
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      character*3 ubi_l
      integer :: ng(3)
      integer :: nv(3)
      ubi_l=adjustl(ubi)
c
#ifdef USE_HDF5
      call fopen(trim(allfile)//'.rst',fid_rst)
      nv = 0
      ng(1) = size(ro_core,1)
      ng(2) = size(ro_core,2)
      ng(3) = size(ro_core,3)
      call fread(MPI_COMM_WORLD,ro_core,ng,nv,"ro_core/"//
     &           trim(ubi_l),fid_rst)
      call fclose(fid_rst)
#else
      if (maswrk) then
        open(3,file=trim(allfile)//'_ro_core_'//trim(ubi_l)//'.rst',
     &       status='old',form="unformatted")
        read(3)ro_core
        close(3)
      endif
c
#ifdef MPI
      if (goparr) then
        call brdcst(msgdbl,ro_core,8*(nrmax+1)*nsort*nspin_0,master,
     &              MPI_COMM_WORLD)
      endif
#endif
#endif
      end
c
c-----------------------------------------------------------------------
c
      subroutine read_wfn_effhf_rst
c
c     Read the wave function and effective Hamiltonian data
c
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use files_mod
#ifdef USE_HDF5
      use hdf5io_mod
#endif
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      integer :: mpierr
      integer :: i_len,ind_k,k0,ispin,iatom,ii
      logical :: oexists
      character*31 :: k0_ch
      character*31 :: me3_k_ch
      character*31 :: me3_tau_ch
      character*3 ubi_l
      integer :: ng(7)
      integer :: nv(7)
      ubi_l=adjustl(ubi)

      i_len=len_trim(allfile)
#ifdef USE_HDF5
      call fopen(trim(allfile)//'.rst',fid_rst)
      if (ubi.eq.'dft'.and.iter_dft<iter_max_scf.and.
     $    iter.eq.iter_dft+1) then
        if(.not.fexist(fid_rst,allfile(1:i_len)//'/v_h_bb/'//
     $                trim(ubi_l))) then
c
c         Apparently the available data is from a DFT calculation
c         but some of the data needed for subsequent calculations
c         such as Hartree-Fock, GW or QP were not produced. Hence
c         repeat the last DFT iteration to construct this data.
c
          iter = iter_dft
        endif
      endif
      nv = 0
      ng(1) = size(eny,1)
      ng(2) = size(eny,2)
      ng(3) = size(eny,3)
      ng(4) = size(eny,4)
      call fread(MPI_COMM_WORLD,eny,ng,nv,
     &           allfile(1:i_len)//'/eny/'//trim(ubi_l),fid_rst)
      ng(1) = size(e_bnd,1)
      ng(2) = size(e_bnd,2)
      ng(3) = size(e_bnd,3)
      call fread(MPI_COMM_WORLD,e_bnd,ng,nv,
     &           allfile(1:i_len)//'/e_bnd/'//trim(ubi_l),fid_rst)
      ng(1) = size(n_bnd,1)
      ng(2) = size(n_bnd,2)
      call fread(MPI_COMM_WORLD,n_bnd,ng,nv,
     &           allfile(1:i_len)//'/n_bnd/'//trim(ubi_l),fid_rst)
      ng(1) = size(e_core,1)
      ng(2) = size(e_core,2)
      ng(3) = size(e_core,3)
      call fread(MPI_COMM_WORLD,e_core,ng,nv,
     &           allfile(1:i_len)//'/e_core/'//trim(ubi_l),fid_rst)
      ng(1) = size(cor_norm,1)
      ng(2) = size(cor_norm,2)
      ng(3) = size(cor_norm,3)
      call fread(MPI_COMM_WORLD,cor_norm,ng,nv,
     &           allfile(1:i_len)//'/cor_norm/'//trim(ubi_l),fid_rst)
      ng(1) = size(pcor,1)
      call fread(MPI_COMM_WORLD,pcor,ng(1),nv(1),
     &           allfile(1:i_len)//'/pcor/'//trim(ubi_l),fid_rst)
      ng(1) = size(pcor_old,1)
      call fread(MPI_COMM_WORLD,pcor_old,ng(1),nv(1),
     &           allfile(1:i_len)//'/pcor_old/'//trim(ubi_l),fid_rst)
      ng(1) = size(qcor,1)
      call fread(MPI_COMM_WORLD,qcor,ng(1),nv(1),
     &           allfile(1:i_len)//'/qcor/'//trim(ubi_l),fid_rst)
      ng(1) = size(qcor_old,1)
      call fread(MPI_COMM_WORLD,qcor_old,ng(1),nv(1),
     &           allfile(1:i_len)//'/qcor_old/'//trim(ubi_l),fid_rst)
      ng(1) = size(p_f,1)
      ng(2) = size(p_f,2)
      ng(3) = size(p_f,3)
      ng(4) = size(p_f,4)
      ng(5) = size(p_f,5)
      call fread(MPI_COMM_WORLD,p_f,ng,nv,
     &           allfile(1:i_len)//'/p_f/'//trim(ubi_l),fid_rst)
      ng(1) = size(q_f,1)
      ng(2) = size(q_f,2)
      ng(3) = size(q_f,3)
      ng(4) = size(q_f,4)
      ng(5) = size(q_f,5)
      call fread(MPI_COMM_WORLD,q_f,ng,nv,
     &           allfile(1:i_len)//'/q_f/'//trim(ubi_l),fid_rst)
      ng(1) = size(pd_f,1)
      ng(2) = size(pd_f,2)
      ng(3) = size(pd_f,3)
      ng(4) = size(pd_f,4)
      ng(5) = size(pd_f,5)
      call fread(MPI_COMM_WORLD,pd_f,ng,nv,
     &           allfile(1:i_len)//'/pd_f/'//trim(ubi_l),fid_rst)
      ng(1) = size(qd_f,1)
      ng(2) = size(qd_f,2)
      ng(3) = size(qd_f,3)
      ng(4) = size(qd_f,4)
      ng(5) = size(qd_f,5)
      call fread(MPI_COMM_WORLD,qd_f,ng,nv,
     &           allfile(1:i_len)//'/qd_f/'//trim(ubi_l),fid_rst)
      ng(1) = size(pd2_f,1)
      ng(2) = size(pd2_f,2)
      ng(3) = size(pd2_f,3)
      ng(4) = size(pd2_f,4)
      ng(5) = size(pd2_f,5)
      call fread(MPI_COMM_WORLD,pd2_f,ng,nv,
     &           allfile(1:i_len)//'/pd2_f/'//trim(ubi_l),fid_rst)
      ng(1) = size(qd2_f,1)
      ng(2) = size(qd2_f,2)
      ng(3) = size(qd2_f,3)
      ng(4) = size(qd2_f,4)
      ng(5) = size(qd2_f,5)
      call fread(MPI_COMM_WORLD,qd2_f,ng,nv,
     &           allfile(1:i_len)//'/qd2_f/'//trim(ubi_l),fid_rst)
      ng(1) = size(gfun,1)
      ng(2) = size(gfun,2)
      call fread(MPI_COMM_WORLD,gfun,ng,nv,
     &           allfile(1:i_len)//'/gfun/'//trim(ubi_l),fid_rst)
      ng(1) = size(gfund,1)
      ng(2) = size(gfund,2)
      call fread(MPI_COMM_WORLD,gfund,ng,nv,
     &           allfile(1:i_len)//'/gfund/'//trim(ubi_l),fid_rst)
      ng(1) = size(gfun_old,1)
      ng(2) = size(gfun_old,2)
      call fread(MPI_COMM_WORLD,gfun_old,ng,nv,
     &           allfile(1:i_len)//'/gfun_old/'//trim(ubi_l),fid_rst)
      ng(1) = size(gfund_old,1)
      ng(2) = size(gfund_old,2)
      call fread(MPI_COMM_WORLD,gfund_old,ng,nv,
     &           allfile(1:i_len)//'/gfund_old/'//trim(ubi_l),fid_rst)
      ng(1) = size(fi,1)
      ng(2) = size(fi,2)
      ng(3) = size(fi,3)
      ng(4) = size(fi,4)
      call fread(MPI_COMM_WORLD,fi,ng,nv,
     &           allfile(1:i_len)//'/fi/'//trim(ubi_l),fid_rst)
      ng(1) = size(dfi,1)
      ng(2) = size(dfi,2)
      ng(3) = size(dfi,3)
      ng(4) = size(dfi,4)
      call fread(MPI_COMM_WORLD,dfi,ng,nv,
     &           allfile(1:i_len)//'/dfi/'//trim(ubi_l),fid_rst)
      ng(1) = size(fidot,1)
      ng(2) = size(fidot,2)
      ng(3) = size(fidot,3)
      ng(4) = size(fidot,4)
      call fread(MPI_COMM_WORLD,fidot,ng,nv,
     &           allfile(1:i_len)//'/fidot/'//trim(ubi_l),fid_rst)
      ng(1) = size(dfidot,1)
      ng(2) = size(dfidot,2)
      ng(3) = size(dfidot,3)
      ng(4) = size(dfidot,4)
      call fread(MPI_COMM_WORLD,dfidot,ng,nv,
     &           allfile(1:i_len)//'/dfidot/'//trim(ubi_l),fid_rst)
      ng(1) = size(fidot2,1)
      ng(2) = size(fidot2,2)
      ng(3) = size(fidot2,3)
      ng(4) = size(fidot2,4)
      call fread(MPI_COMM_WORLD,fidot2,ng,nv,
     &           allfile(1:i_len)//'/fidot2/'//trim(ubi_l),fid_rst)
      ng(1) = size(dfidot2,1)
      ng(2) = size(dfidot2,2)
      ng(3) = size(dfidot2,3)
      ng(4) = size(dfidot2,4)
      call fread(MPI_COMM_WORLD,dfidot2,ng,nv,
     &           allfile(1:i_len)//'/dfidot2/'//trim(ubi_l),fid_rst)
      ng(1) = size(ptnl,1)
      ng(2) = size(ptnl,2)
      ng(3) = size(ptnl,3)
      ng(4) = size(ptnl,4)
      call fread(MPI_COMM_WORLD,ptnl,ng,nv,
     &           allfile(1:i_len)//'/ptnl/'//trim(ubi_l),fid_rst)
      ng(1) = size(ffsmt,1)
      ng(2) = size(ffsmt,2)
      ng(3) = size(ffsmt,3)
      ng(4) = size(ffsmt,4)
      ng(5) = size(ffsmt,5)
      ng(6) = size(ffsmt,6)
      ng(7) = size(ffsmt,7)
      call fread(MPI_COMM_WORLD,ffsmt,ng,nv,
     &           allfile(1:i_len)//'/ffsmt/'//trim(ubi_l),fid_rst)
      ng(1) = size(v_mt_h,1)
      call fread(MPI_COMM_WORLD,v_mt_h,ng(1),nv(1),
     &           allfile(1:i_len)//'/v_mt_h/'//trim(ubi_l),fid_rst)
      ng(1) = size(v_mt_xc,1)
      call fread(MPI_COMM_WORLD,v_mt_xc,ng(1),nv(1),
     &           allfile(1:i_len)//'/v_mt_xc/'//trim(ubi_l),fid_rst)
      ng(1) = size(v_intr_h,1)
      call fread(MPI_COMM_WORLD,v_intr_h,ng(1),nv(1),
     &           allfile(1:i_len)//'/v_intr_h/'//trim(ubi_l),fid_rst)
      ng(1) = size(v_intr_xc,1)
      ng(2) = size(v_intr_xc,2)
      call fread(MPI_COMM_WORLD,v_intr_xc,ng,nv,
     &           allfile(1:i_len)//'/v_intr_xc/'//trim(ubi_l),fid_rst)
      if (ubi.eq.' gw'.or.ubi.eq.' qp') then
        ng(1) = size(e_qp,1)
        ng(2) = size(e_qp,2)
        ng(3) = size(e_qp,3)
        call fread(MPI_COMM_WORLD,e_qp,ng,nv,
     &             allfile(1:i_len)//'/e_qp/'//trim(ubi_l),fid_rst)
      endif
c
c     Distributed data structures
c     - distributed over Tau
c     - replicated within k
c
      nv(3) = n3_mpi_k(me3_k+1)
      ng(1) = size(ev_bnd,1)
      ng(2) = size(ev_bnd,2)
      ng(3) = npnt
      ng(4) = size(ev_bnd,4)
      do ii = 1, ng(4)
        nv(4) = ii-1
        call fread(comm_k,ev_bnd(:,:,1:ndim3_k(me3_k+1),ii:ii),ng,nv,
     &             allfile(1:i_len)//'/ev_bnd/'//trim(ubi_l),fid_rst)
      enddo
      nv(3) = n3_mpi_k(me3_k+1)
      ng(1) = size(z_bnd,1)
      ng(2) = size(z_bnd,2)
      ng(3) = npnt
      ng(4) = size(z_bnd,4)
      do ii = 1, ng(4)
        nv(4) = ii-1
        call fread(comm_k,z_bnd(:,:,1:ndim3_k(me3_k+1),ii:ii),ng,nv,
     &             allfile(1:i_len)//'/z_bnd/'//trim(ubi_l),fid_rst)
      enddo
      if (ubi.ne.'dft') then
        nv(3) = n3_mpi_k(me3_k+1)
        ng(1) = size(v_h_bb,1)
        ng(2) = size(v_h_bb,2)
        ng(3) = npnt
        ng(4) = size(v_h_bb,4)
        do ii = 1, ng(4)
          nv(4) = ii-1
          call fread(comm_k,v_h_bb(:,:,1:ndim3_k(me3_k+1),ii:ii),ng,nv,
     &               allfile(1:i_len)//'/v_h_bb/'//trim(ubi_l),fid_rst)
        enddo
        nv(3) = n3_mpi_k(me3_k+1)
        ng(1) = size(v_h_bb_old,1)
        ng(2) = size(v_h_bb_old,2)
        ng(3) = npnt
        ng(4) = size(v_h_bb_old,4)
        do ii = 1, ng(4)
          nv(4) = ii-1
          call fread(comm_k,v_h_bb_old(:,:,1:ndim3_k(me3_k+1),ii:ii),
     &               ng,nv,
     &               allfile(1:i_len)//'/v_h_bb_old/'//trim(ubi_l),
     &               fid_rst)
        enddo
        nv(3) = n3_mpi_k(me3_k+1)
        ng(1) = size(v_xc_bb,1)
        ng(2) = size(v_xc_bb,2)
        ng(3) = npnt
        ng(4) = size(v_xc_bb,4)
        do ii = 1, ng(4)
          nv(4) = ii-1
          call fread(comm_k,v_xc_bb(:,:,1:ndim3_k(me3_k+1),ii:ii),ng,nv,
     &               allfile(1:i_len)//'/v_xc_bb/'//trim(ubi_l),fid_rst)
        enddo
        nv(3) = n3_mpi_k(me3_k+1)
        ng(1) = size(v_xc_bb_old,1)
        ng(2) = size(v_xc_bb_old,2)
        ng(3) = npnt
        ng(4) = size(v_xc_bb_old,4)
        do ii = 1, ng(4)
          nv(4) = ii-1
          call fread(comm_k,v_xc_bb_old(:,:,1:ndim3_k(me3_k+1),ii:ii),
     &               ng,nv,
     &               allfile(1:i_len)//'/v_xc_bb_old/'//trim(ubi_l),
     &               fid_rst)
        enddo
        nv(3) = n3_mpi_k(me3_k+1)
        ng(1) = size(sigx_solid_k,1)
        ng(2) = size(sigx_solid_k,2)
        ng(3) = npnt
        ng(4) = size(sigx_solid_k,4)
        do ii = 1, ng(4)
          nv(4) = ii-1
          call fread(comm_k,sigx_solid_k(:,:,1:ndim3_k(me3_k+1),ii:ii),
     &               ng,nv,
     &               allfile(1:i_len)//'/sigx_solid_k/'//trim(ubi_l),
     &               fid_rst)
        enddo
        nv(3) = n3_mpi_k(me3_k+1)
        ng(1) = size(sigx_solid_k_old,1)
        ng(2) = size(sigx_solid_k_old,2)
        ng(3) = npnt
        ng(4) = size(sigx_solid_k_old,4)
        do ii = 1, ng(4)
          nv(4) = ii-1
          call fread(comm_k,
     &               sigx_solid_k_old(:,:,1:ndim3_k(me3_k+1),ii:ii),
     &               ng,nv,
     &               allfile(1:i_len)//'/sigx_solid_k_old/'
     &                      //trim(ubi_l),
     &               fid_rst)
        enddo
      endif
      if (ubi.eq.' gw'.or.ubi.eq.' qp') then
        nv(3) = n3_mpi_k(me3_k+1)
        ng(1) = size(z_qp,1)
        ng(2) = size(z_qp,2)
        ng(3) = npnt
        ng(4) = size(z_qp,4)
        do ii = 1, ng(4)
          nv(4) = ii-1
          call fread(comm_k,z_qp(:,:,1:ndim3_k(me3_k+1),ii:ii),ng,nv,
     &               allfile(1:i_len)//'/z_qp/'//trim(ubi_l),
     &               fid_rst)
        enddo
        nv(3) = n3_mpi_k(me3_k+1)
        ng(1) = size(ev_qp,1)
        ng(2) = size(ev_qp,2)
        ng(3) = npnt
        ng(4) = size(ev_qp,4)
        do ii = 1, ng(4)
          nv(4) = ii-1
          call fread(comm_k,ev_qp(:,:,1:ndim3_k(me3_k+1),ii:ii),ng,nv,
     &               allfile(1:i_len)//'/ev_qp/'//trim(ubi_l),
     &               fid_rst)
        enddo
        nv(3) = n3_mpi_k(me3_k+1)
        ng(1) = size(q_qp,1)
        ng(2) = size(q_qp,2)
        ng(3) = npnt
        ng(4) = size(q_qp,4)
        do ii = 1, ng(4)
          nv(4) = ii-1
          call fread(comm_k,q_qp(:,:,1:ndim3_k(me3_k+1),ii:ii),ng,nv,
     &               allfile(1:i_len)//'/q_qp/'//trim(ubi_l),
     &               fid_rst)
        enddo
        nv = 0
        nv(4) = size(sig_c_tau_old,4)*me3_tau
        nv(5) = n3_mpi_k(me3_k+1)
        ng(1) = size(sig_c_tau_old,1)
        ng(2) = size(sig_c_tau_old,2)
        ng(3) = size(sig_c_tau_old,3)
        ng(4) = size(sig_c_tau_old,4)*nproc_tau
        ng(5) = npnt
        ng(6) = size(sig_c_tau_old,6)
        do ii = 1, ng(6)
          nv(6) = ii-1
          call fread(comm_tau_k,
     &               sig_c_tau_old(:,:,:,:,1:ndim3_k(me3_k+1),ii:ii),
     &               ng,nv,
     &               allfile(1:i_len)//'/sig_c_tau_old/'//trim(ubi_l),
     &               fid_rst)
        enddo
      endif
      call fclose(fid_rst)
#else
      if (maswrk) then
        open(3,file=allfile(1:i_len)//'_eny_'//trim(ubi_l)//'.rst'
     $    ,status='old',form='unformatted')        
        read(3) eny
        close(3)
        open(3,file=allfile(1:i_len)//'_e_bnd_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) e_bnd
        close(3)
        open(3,file=allfile(1:i_len)//'_n_bnd_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) n_bnd
        close(3)
        open(3,file=allfile(1:i_len)//'_e_core_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) e_core
        close(3)
        open(3,file=allfile(1:i_len)//'_cor_norm_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) cor_norm
        close(3)
        open(3,file=allfile(1:i_len)//'_pcor_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) pcor
        close(3)
        open(3,file=allfile(1:i_len)//'_pcor_old_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) pcor_old
        close(3)
        open(3,file=allfile(1:i_len)//'_qcor_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) qcor
        close(3)
        open(3,file=allfile(1:i_len)//'_qcor_old_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) qcor_old
        close(3)
        open(3,file=allfile(1:i_len)//'_p_f_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) p_f
        close(3)
        open(3,file=allfile(1:i_len)//'_q_f_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) q_f
        close(3)
        open(3,file=allfile(1:i_len)//'_pd_f_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) pd_f
        close(3)
        open(3,file=allfile(1:i_len)//'_qd_f_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) qd_f
        close(3)
        open(3,file=allfile(1:i_len)//'_pd2_f_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) pd2_f
        close(3)
        open(3,file=allfile(1:i_len)//'_qd2_f_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) qd2_f
        close(3)
        open(3,file=allfile(1:i_len)//'_gfun_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) gfun
        close(3)
        open(3,file=allfile(1:i_len)//'_gfund_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) gfund
        close(3)
        open(3,file=allfile(1:i_len)//'_gfun_old_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) gfun_old
        close(3)
        open(3,file=allfile(1:i_len)//'_gfund_old_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) gfund_old
        close(3)
        open(3,file=allfile(1:i_len)//'_fi_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) fi
        close(3)
        open(3,file=allfile(1:i_len)//'_dfi_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) dfi
        close(3)
        open(3,file=allfile(1:i_len)//'_fidot_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) fidot
        close(3)
        open(3,file=allfile(1:i_len)//'_dfidot_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) dfidot
        close(3)
        open(3,file=allfile(1:i_len)//'_fidot2_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) fidot2
        close(3)
        open(3,file=allfile(1:i_len)//'_dfidot2_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) dfidot2
        close(3)
        open(3,file=allfile(1:i_len)//'_ptnl_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) ptnl
        close(3)
        open(3,file=allfile(1:i_len)//'_ffsmt_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) ffsmt
        close(3)
        open(3,file=allfile(1:i_len)//'_v_mt_h_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) v_mt_h
        close(3)
        open(3,file=allfile(1:i_len)//'_v_mt_xc_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) v_mt_xc
        close(3)
        open(3,file=allfile(1:i_len)//'_v_intr_h_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) v_intr_h
        close(3)
        open(3,file=allfile(1:i_len)//'_v_intr_xc_'
     $    //trim(ubi_l)//'.rst',status='old',form='unformatted')        
        read(3) v_intr_xc
        close(3)

        if(ubi.eq.' gw'.or.ubi.eq.' qp') then        
          open(3,file=allfile(1:i_len)//'_e_qp_'
     $      //trim(ubi_l)//'.rst',status='old',form='unformatted')        
          read(3) e_qp
          close(3)
        endif
      endif
#ifdef MPI
      if (goparr) then
        call mpi_bcast(eny,maxntle*(nrel*maxb+1)*nsort*nspin_0,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(e_bnd,nbndf*npnt*nspin_0,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(n_bnd,npnt*nspin_0,
     $    mpi_integer,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(e_core,ncormax*nspin_0*nsort,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(cor_norm,ncormax*nspin_0*nsort,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(pcor,maxmtcor,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(pcor_old,maxmtcor,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(qcor,maxmtcor,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(qcor_old,maxmtcor,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(p_f,
     $    (1+maxnrad)*maxntle*(nrel*maxb+1)*nsort*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(q_f,
     $    (1+maxnrad)*maxntle*(nrel*maxb+1)*nsort*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(pd_f,
     $    (1+maxnrad)*maxntle*(nrel*maxb+1)*nsort*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(qd_f,
     $    (1+maxnrad)*maxntle*(nrel*maxb+1)*nsort*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(pd2_f,
     $    (1+maxnrad)*maxntle*(nrel*maxb+1)*nsort*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(qd2_f,
     $    (1+maxnrad)*maxntle*(nrel*maxb+1)*nsort*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(gfun,
     $    maxwf*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(gfund,
     $    maxwf*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(gfun_old,
     $    maxwf*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(gfund_old,
     $    maxwf*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(fi,
     $    maxntle*(nrel*maxb+1)*nsort*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(dfi,
     $    maxntle*(nrel*maxb+1)*nsort*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(fidot,
     $    maxntle*(nrel*maxb+1)*nsort*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(dfidot,
     $    maxntle*(nrel*maxb+1)*nsort*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(fidot2,
     $    maxntle*(nrel*maxb+1)*nsort*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(dfidot2,
     $    maxntle*(nrel*maxb+1)*nsort*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(ptnl,
     $    maxntle*(maxb+1)*nsort*2,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(ffsmt,
     $    4*maxntle**2*(nrel*maxb+1)*nsort*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)
        
        call mpi_bcast(v_mt_h,maxmt/nspin_0,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(v_mt_xc,maxmt,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(v_intr_h,nplwro,
     $    mpi_double_complex,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(v_intr_xc,nplwro*nspin_0,
     $    mpi_double_complex,0,MPI_COMM_WORLD,mpierr)

        if(ubi.eq.' gw'.or.ubi.eq.' qp') then
          call mpi_bcast(e_qp,nbndf*npnt*nspin_0,
     $      mpi_double_precision,0,MPI_COMM_WORLD,mpierr)
        endif

      endif
#endif     
      if (maswrk) then 
        k0_ch = "1"
        k0_ch = adjustl(k0_ch)
        if (ubi.eq.'dft'.and.iter_dft<iter_max_scf.and.
     $      iter.eq.iter_dft+1) then
          inquire(file=allfile(1:i_len)//'_v_h_bb_'//trim(k0_ch)
     $      //'_'//trim(ubi_l)//'.rst',exist=oexists)
          if (.not.oexists) then
c
c           Apparently the available data is from a DFT calculation
c           but some of the data needed for subsequent calculations
c           such as Hartree-Fock, GW or QP were not produced. Hence
c           repeat the last DFT iteration to construct this data.
c
            iter = iter_dft
          endif
        endif
      endif
#ifdef MPI
      if (goparr) then
        call mpi_bcast(iter,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      endif
#endif
      if (me3_tau .eq. 0) then
        do ind_k=1,ndim3_k(me3_k+1)
          k0=n3_mpi_k(me3_k+1)+ind_k
          write(k0_ch,*) k0
          k0_ch = adjustl(k0_ch)
          
          open(3,file=allfile(1:i_len)//'_ev_bnd_'//trim(k0_ch)
     $      //'_'//trim(ubi_l)//'.rst',status='old',form='unformatted')
          do ispin=1, nspin_0
            read(3) ev_bnd(:,:,ind_k,ispin)
          enddo
          close(3)
          open(3,file=allfile(1:i_len)//'_z_bnd_'//trim(k0_ch)
     $      //'_'//trim(ubi_l)//'.rst',status='old',form='unformatted')
          do ispin=1, nspin_0
            read(3) z_bnd(:,:,ind_k,ispin)
          enddo
          close(3)

          if (ubi.ne.'dft') then
            open(3,file=allfile(1:i_len)//'_v_h_bb_'//trim(k0_ch)
     $        //'_'//trim(ubi_l)//'.rst',status='old',
     $        form='unformatted')        
            do ispin=1, nspin_0
              read(3) v_h_bb(:,:,ind_k,ispin)
            enddo
            close(3)
            open(3,file=allfile(1:i_len)//'_v_h_bb_old_'//trim(k0_ch)
     $        //'_'//trim(ubi_l)//'.rst',status='old',
     $        form='unformatted')        
            do ispin=1, nspin_0
              read(3) v_h_bb_old(:,:,ind_k,ispin)
            enddo
            close(3)
            open(3,file=allfile(1:i_len)//'_v_xc_bb_'//trim(k0_ch)
     $        //'_'//trim(ubi_l)//'.rst',status='old',
     $        form='unformatted')        
            do ispin=1, nspin_0
              read(3) v_xc_bb(:,:,ind_k,ispin)
            enddo
            close(3)
            open(3,file=allfile(1:i_len)//'_v_xc_bb_old_'//trim(k0_ch)
     $        //'_'//trim(ubi_l)//'.rst',status='old',
     $        form='unformatted')        
            do ispin=1, nspin_0
              read(3) v_xc_bb_old(:,:,ind_k,ispin)
            enddo
            close(3)

            open(3,file=allfile(1:i_len)//'_sigx_solid_k_'//
     $        trim(k0_ch)//'_'//trim(ubi_l)//'.rst',status='old',
     $        form='unformatted')        
            do ispin=1, nspin_0
              read(3) sigx_solid_k(:,:,ind_k,ispin)
            enddo
            close(3)
            open(3,file=allfile(1:i_len)//'_sigx_solid_k_old_'//
     $        trim(k0_ch)//'_'//trim(ubi_l)//'.rst',status='old',
     $        form='unformatted')        
            do ispin=1, nspin_0
              read(3) sigx_solid_k_old(:,:,ind_k,ispin)
            enddo
            close(3)
            if(ubi.eq.' gw'.or.ubi.eq.' qp') then
              open(3,file=allfile(1:i_len)//'_z_qp_'//trim(k0_ch)
     $          //'_'//trim(ubi_l)//'.rst',status='old',
     $          form='unformatted')        
              do ispin=1, nspin_0
                read(3) z_qp(:,:,ind_k,ispin)
              enddo
              close(3)
              open(3,file=allfile(1:i_len)//'_ev_qp_'//trim(k0_ch)
     $          //'_'//trim(ubi_l)//'.rst',status='old',
     $          form='unformatted')        
              do ispin=1, nspin_0
                read(3) ev_qp(:,:,ind_k,ispin)
              enddo
              close(3)
              open(3,file=allfile(1:i_len)//'_q_qp_'//trim(k0_ch)
     $          //'_'//trim(ubi_l)//'.rst',status='old',
     $          form='unformatted')        
              do ispin=1, nspin_0
                read(3) q_qp(:,:,ind_k,ispin)
              enddo
              close(3)
            endif
          endif
        enddo
      endif
      if (ubi .eq. ' gw'.or. ubi .eq. ' qp') then
        write(me3_k_ch,*)me3_k
        me3_k_ch = adjustl(me3_k_ch)
        write(me3_tau_ch,*)me3_tau
        me3_tau_ch = adjustl(me3_tau_ch)
        open(3,file=allfile(1:i_len)//'_sig_c_tau_old_'//trim(me3_k_ch)
     $    //'_'//trim(me3_tau_ch)//'_'//trim(ubi_l)//'.rst',
     $    status='old',form='unformatted')        
        read(3) sig_c_tau_old
        close(3)
      endif

#ifdef MPI
      if (goparr) then
c
        call mpi_bcast(ev_bnd,nbasmpw*nbndf*maxval(ndim3_k)*nspin_0,
     $    mpi_double_complex,0,comm_k,mpierr)

        call mpi_bcast(z_bnd,nfun*nbndf*maxval(ndim3_k)*nspin_0,
     $    mpi_double_complex,0,comm_k,mpierr)

        if (ubi.ne.'dft') then
          call mpi_bcast(v_h_bb_old,
     $      nbndf*nbndf*maxval(ndim3_k)*nspin_0,
     $      mpi_double_precision,0,comm_k,mpierr)
          call mpi_bcast(v_h_bb,nbndf*nbndf*maxval(ndim3_k)*nspin_0,
     $      mpi_double_precision,0,comm_k,mpierr)
          call mpi_bcast(v_xc_bb_old,
     $      nbndf*nbndf*maxval(ndim3_k)*nspin_0,
     $      mpi_double_precision,0,comm_k,mpierr)
          call mpi_bcast(v_xc_bb,nbndf*nbndf*maxval(ndim3_k)*nspin_0,
     $      mpi_double_precision,0,comm_k,mpierr)
          call mpi_bcast(sigx_solid_k_old,
     $      nbndf*nbndf*maxval(ndim3_k)*nspin_0,
     $      mpi_double_precision,0,comm_k,mpierr)
          call mpi_bcast(sigx_solid_k,
     $      nbndf*nbndf*maxval(ndim3_k)*nspin_0,
     $      mpi_double_precision,0,comm_k,mpierr)
          if(ubi.eq.' gw'.or.ubi.eq.' qp') then
            call mpi_bcast(z_qp,nfun*nbndf*maxval(ndim3_k)*nspin_0,
     $        mpi_double_complex,0,comm_k,mpierr)            
            call mpi_bcast(ev_qp,nbasmpw*nbndf*maxval(ndim3_k)*
     $        nspin_0,mpi_double_complex,0,comm_k,mpierr)
            call mpi_bcast(q_qp,nbndf*nbndf*maxval(ndim3_k)*nspin_0,
     $        mpi_double_complex,0,comm_k,mpierr)
          endif
        endif
      endif
#endif
#endif
      
      if (maswrk) write(iun,'(a, a,a)')
     $  'reading wfn from restart file in ', ubi, ' step'
      end
c
c-----------------------------------------------------------------------
c
      subroutine read_green_effhf_rst
c
c     Read the green function from the restart files
c
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use files_mod
#ifdef USE_HDF5
      use hdf5io_mod
#endif
c
      implicit none
c
#ifdef MPI
      include 'mpif.h'
#endif
c
      integer :: i_len,ind_k,k0,ispin
      integer :: ii
      integer :: mpierr
      character*31 :: k0_ch
      character*31 :: me3_k_ch
      character*31 :: me3_tau_ch
      character*3  :: ubi_l
      integer :: ng(6)
      integer :: nv(6)
      logical :: file_exists
c
      ubi_l=adjustl(ubi)
      i_len=len_trim(allfile)
#ifdef USE_HDF5
      call fopen(trim(allfile)//'.rst',fid_rst)
      nv = 0
      ng(1) = size(g_loc_0,1)
      ng(2) = size(g_loc_0,2)
      ng(3) = size(g_loc_0,3)
      ng(4) = size(g_loc_0,4)
      call fread(MPI_COMM_WORLD,g_loc_0,ng,nv,
     &           allfile(1:i_len)//'/g_loc_0/'//trim(ubi_l),fid_rst)
      if ((ubi .eq. ' gw') .or. (ubi .eq. ' qp')) then
        nv(3) = n3_mpi_k(me3_k+1)
        ng(1) = size(z_ren,1)
        ng(2) = size(z_ren,2)
        ng(3) = npnt
        ng(4) = size(z_ren,4)
        do ii = 1, ng(4)
          nv(4) = ii-1
          call fread(comm_k,z_ren(:,:,1:ndim3_k(me3_k+1),ii:ii),ng,nv,
     &               allfile(1:i_len)//'/z_ren/'//trim(ubi_l),fid_rst)
        enddo
      endif
      if ((ubi .eq. 'dft') .or. (ubi .eq. ' hf')) then
        nv    = 0
        nv(2) = n3_mpi_k(me3_k+1)
        ng(1) = size(g_full_00,1)
        ng(2) = npnt
        ng(3) = size(g_full_00,3)
        do ii = 1, ng(3)
          nv(3) = ii-1
          call fread(comm_k,g_full_00(:,1:ndim3_k(me3_k+1),ii:ii),ng,nv,
     &               allfile(1:i_len)//'/g_full_00/'//trim(ubi_l),
     &               fid_rst)
        enddo
      endif
      if (ubi .eq. ' gw' .or. (ubi .eq. ' qp')) then
        nv    = 0
        nv(3) = n3_mpi_k(me3_k+1)
        ng(1) = size(g_full_0,1)
        ng(2) = size(g_full_0,2)
        ng(3) = npnt
        ng(4) = size(g_full_0,4)
        do ii = 1, ng(4)
          nv(4) = ii-1
          call fread(comm_k,g_full_0(:,:,1:ndim3_k(me3_k+1),ii:ii),
     &               ng,nv,
     &               allfile(1:i_len)//'/g_full_0/'//trim(ubi_l),
     &               fid_rst)
        enddo
c
        nv    = 0
        nv(4) = size(g_full,4)*me3_tau
        nv(5) = n3_mpi_k(me3_k+1)
        ng(1) = size(g_full,1)
        ng(2) = size(g_full,2)
        ng(3) = size(g_full,3)
        ng(4) = size(g_full,4)*nproc_tau
        ng(5) = npnt
        ng(6) = size(g_full,6)
        do ii = 1, ng(6)
          nv(6) = ii-1
          call fread(comm_tau_k,
     &               g_full(:,:,:,:,1:ndim3_k(me3_k+1),ii:ii),ng,nv,
     &               allfile(1:i_len)//'/g_full/'//trim(ubi_l),
     &               fid_rst)
        enddo
      endif
      call fclose(fid_rst)
#else
      if (maswrk) then
          open(3,file=allfile(1:i_len)//'_g_loc_0_'
     $      //trim(ubi_l)//'.rst',status='old',form='unformatted')
          read(3) g_loc_0
          close(3)
      endif
c
#ifdef MPI
      if (goparr) then
          call mpi_bcast(g_loc_0,maxel**2*natom*nspin_0,
     $                   mpi_double_complex,0,MPI_COMM_WORLD,mpierr)
      endif
#endif
c
      if (me3_tau .eq. 0) then
        do ind_k=1,ndim3_k(me3_k+1)
          k0=n3_mpi_k(me3_k+1)+ind_k
          write(k0_ch,*)k0
          k0_ch=adjustl(k0_ch)

          if ((ubi .eq. ' gw') .or. (ubi .eq. ' qp')) then
            open(3,file=allfile(1:i_len)//'_z_ren_'//trim(k0_ch)
     $        //'_'//trim(ubi_l)//'.rst',status='old',
     $        form='unformatted')        
            do ispin=1, nspin_0
              read(3) z_ren(:,:,ind_k,ispin)
            enddo
            close(3)
          endif
          
          if ((ubi .eq. 'dft') .or. (ubi .eq. ' hf')) then
            open(3,file=allfile(1:i_len)//'_g_full_00_'//trim(k0_ch)
     $        //'_'//trim(ubi_l)//'.rst',status='old',
     $        form='unformatted')        
            do ispin=1, nspin_0
              read(3) g_full_00(:,ind_k, ispin)
            enddo
            close(3)
          endif

          if (ubi_0 .eq. ' gw' .or. (ubi .eq. ' qp')) then
            open(3,file=allfile(1:i_len)//'_g_full_0_'//trim(k0_ch)
     $        //'_'//trim(ubi_l)//'.rst',status='old',
     $        form='unformatted')        
            do ispin=1, nspin_0
              read(3) g_full_0(:,:,ind_k, ispin)
            enddo
            close(3)
          endif
        enddo
      endif
#ifdef MPI
      if (goparr) then
        if  ((ubi .eq. ' gw') .or. (ubi .eq. ' qp')) then
          call mpi_bcast(z_ren,nbndf**2*maxval(ndim3_k)*nspin_0,
     $      mpi_double_precision,0,comm_k,mpierr)
        endif
        
        if ((ubi .eq. 'dft') .or. (ubi .eq. ' hf')) then
          call mpi_bcast(g_full_00,nbndf*maxval(ndim3_k)*nspin_0,
     $      mpi_double_precision,0,comm_k,mpierr)
        endif
        if (ubi .eq. ' gw' .or. (ubi .eq. ' qp')) then
          call mpi_bcast(g_full_0,nbndf**2*maxval(ndim3_k)*nspin_0,
     $      mpi_double_complex,0,comm_k,mpierr)
        endif
      endif
#endif
      if (ubi .eq. ' gw' .or. (ubi .eq. ' qp')) then
        write(me3_k_ch,*)me3_k
        me3_k_ch = adjustl(me3_k_ch)
        write(me3_tau_ch,*)me3_tau
        me3_tau_ch = adjustl(me3_tau_ch)
        open(3,file=allfile(1:i_len)//'_g_full_tau_'//trim(me3_k_ch)
     $    //'_'//trim(me3_tau_ch)//'_'//trim(ubi_l)//'.rst',
     $    status='old',form='unformatted')        
        read(3) g_full
        close(3)
      endif
#endif

      if (maswrk) write(iun,'(a, a,a)')
     $  'reading green from restart file in ', ubi, ' step'
      end
c
c-----------------------------------------------------------------------

