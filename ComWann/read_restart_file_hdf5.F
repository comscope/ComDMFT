!     restart
!     1) dft case: read ro at iter_dft=2
!     2) hf case: a_store_ubn is not 
!     (1) read e_bnd, z_bnd, and ev_bnd
!     (1) read v_h_bb_old and v_h_bb, v_xc_bb, v_xc_bb_old
!     (1) read sigx_solid_k_old and sigx_solid_k
!     A_store_ubn: should be selected at each iteration !!! revise it!!!


      subroutine read_chemical_potential_rst_hdf5
      use comwann_mod
      use hdf5io_mod
      use hdf5
      implicit none
      include 'mpif.h'
      integer :: i_len,iread,rst_len
      integer(hid_t) :: fid_rst = -1
      rst_len=len_trim(rstpath)
      call fopen(rstpath(1:rst_len)//'/'//
     $  trim(allfile)//'.rst',fid_rst) 
      call fread(MPI_COMM_WORLD,chem_pot,"chemical_potential/"//
     &           trim(ubi_0)//"/chem_pot",fid_rst)
      if (ubi_0 .eq. 'qp') then
        call fread(MPI_COMM_WORLD,chem_pot_qp,"chemical_potential/"//
     &             trim(ubi_0)//"/chem_pot_qp",fid_rst)
      endif   
      call fclose(fid_rst)   

      if (maswrk) then 
        write(993,'(a, a,a)')
     $  'reading chemical potential from restart file in ',
     $  ubi_0, ' step'
        write(993,*)'chem_pot',chem_pot
        write(iun,*) 'chem_pot_qp', chem_pot_qp
      endif
      end


      subroutine read_wfn_effhf_rst_hdf5
************************************************************
*     store charge density at given iteration              *
************************************************************
      use comwann_mod
      use hdf5io_mod
      use hdf5
      implicit none
      include 'mpif.h'
      integer :: i_len,ind_k,k0,ispin,iatom,ii,rst_len
      character*8 :: k0_ch
      integer(hid_t) :: fid_rst = -1
      integer :: ng(7)
      integer :: nv(7)
      
      rst_len=len_trim(rstpath)
      i_len=len_trim(allfile)
      call fopen(rstpath(1:rst_len)//'/'//
     $ trim(allfile)//'.rst',fid_rst)
      nv    = 0
      ng(1) = size(e_bnd,1)
      ng(2) = size(e_bnd,2)
      ng(3) = size(e_bnd,3)
      call fread(MPI_COMM_WORLD,e_bnd,ng,nv,
     &           allfile(1:i_len)//'/e_bnd/'//trim(ubi_0),fid_rst)
!      write(993,*)'e_bnd',e_bnd
c      if (maswrk) write(180,*)'e_bnd',e_bnd
      ng(1) = size(n_bnd,1)
      ng(2) = size(n_bnd,2)
      call fread(MPI_COMM_WORLD,n_bnd,ng,nv,
     &           allfile(1:i_len)//'/n_bnd/'//trim(ubi_0),fid_rst)
!      if(maswrk)  write(180,*)'n_bnd',n_bnd

      ng(1) = size(gfun,1)
      ng(2) = size(gfun,2)
      call fread(MPI_COMM_WORLD,gfun,ng,nv,
     &           allfile(1:i_len)//'/gfun/'//trim(ubi_0),fid_rst)
!      if(maswrk) write(180,*)'gfun',gfun

      ng(1) = size(gfund,1)
      ng(2) = size(gfund,2)
      call fread(MPI_COMM_WORLD,gfund,ng,nv,
     &           allfile(1:i_len)//'/gfund/'//trim(ubi_0),fid_rst)
!      if(maswrk) write(180,*)'gfund',gfund

      ng(1) = size(ffsmt,1)
      ng(2) = size(ffsmt,2)
      ng(3) = size(ffsmt,3)
      ng(4) = size(ffsmt,4)
      ng(5) = size(ffsmt,5)
      ng(6) = size(ffsmt,6)
      ng(7) = size(ffsmt,7)
      call fread(MPI_COMM_WORLD,ffsmt,ng,nv,
     &           allfile(1:i_len)//'/ffsmt/'//trim(ubi_0),fid_rst)
!      if (maswrk) write(180,*)'ffsmt',ffsmt

      if (ubi_0 .eq. 'qp') then
!        if(maswrk) write(180,*)'ubi_0',ubi_0
        ng(1) = size(e_qp,1)
        ng(2) = size(e_qp,2)
        ng(3) = size(e_qp,3)
        call fread(MPI_COMM_WORLD,e_qp,ng,nv,
     &             allfile(1:i_len)//'/e_qp/'//trim(ubi_0),fid_rst)
      endif

!       if (maswrk) write(180,*)'e_qp',e_qp

      nv(3) = n_mpi_k(me+1)
      ng(1) = size(ev_bnd,1)
      ng(2) = size(ev_bnd,2)
      ng(3) = npnt
      ng(4) = size(ev_bnd,4)
      call fread(MPI_COMM_WORLD,ev_bnd,ng,nv,
     &           allfile(1:i_len)//'/ev_bnd/'//trim(ubi_0),fid_rst)
      if(maswrk) write(iun,*)
     & allfile(1:i_len)//'/ev_bnd/'//trim(ubi_0)
      if(maswrk) write(iun,*)'ev_bnd',ev_bnd
!      write(180,*)'ev_bnd',ev_bnd

      nv(3) = n_mpi_k(me+1)
      ng(1) = size(z_bnd,1)
      ng(2) = size(z_bnd,2)
      ng(3) = npnt
      ng(4) = size(z_bnd,4)
      call fread(MPI_COMM_WORLD,z_bnd,ng,nv,
     &           allfile(1:i_len)//'/z_bnd/'//trim(ubi_0),fid_rst)
!      write(180,*)'z_bnd',z_bnd

      if (ubi_0 .eq. 'qp') then
        nv(3) = n_mpi_k(me+1)
        ng(1) = size(z_qp,1)
        ng(2) = size(z_qp,2)
        ng(3) = npnt
        ng(4) = size(z_qp,4)
        call fread(MPI_COMM_WORLD,z_qp,ng,nv,
     &             allfile(1:i_len)//'/z_qp/'//trim(ubi_0),
     &             fid_rst)
        nv(3) = n_mpi_k(me+1)
        ng(1) = size(ev_qp,1)
        ng(2) = size(ev_qp,2)
        ng(3) = npnt
        ng(4) = size(ev_qp,4)
        call fread(MPI_COMM_WORLD,ev_qp,ng,nv,
     &             allfile(1:i_len)//'/ev_qp/'//trim(ubi_0),
     &             fid_rst)
      endif
!      write(993,*)'ubi_0',ubi_0
chvd  endif
      call mpi_barrier(mpi_comm_world,mpierr)

      call fclose(fid_rst)
      if (maswrk) write(993,'(a, a,a)')
     $  'reading wfn from restart file in ', ubi_0, ' step'
c      write(993,*)'read wfn from restart'
      end

