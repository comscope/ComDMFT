      Program QFT
      use manager_mod
      use parallel_mod
      use solid_mod
      use atom_mod
      use units_mod
      use wannier90_mod
      use vertex_mod
!      use gitinfo_mod

************************************************************
*     The main module                             *
************************************************************

      implicit none
      include 'mpif.h'
      character*300 :: ch
      character*300 :: chbk
      logical :: file_exist
      integer :: i,ispin,k,ind_k,na,nb,nc,ic,ia,ib
      real*8 :: v(3),pi2,fas
      call beging
!      dmft_gw=' QP'
!      allocate(mt2ciw(13,7,2))
!      mt2ciw=0
!      mt2ciw(1,2,1)=1
!      mt2ciw(2,2,1)=2
!      mt2ciw(3,2,1)=3
!      mt2ciw(1,1,2)=4
!      mt2ciw(1,2,2)=5
!      mt2ciw(2,2,2)=6
!      mt2ciw(3,2,2)=7
!      mt2ciw(1,3,2)=8
!      mt2ciw(2,3,2)=9
!      mt2ciw(3,3,2)=10
!      mt2ciw(4,3,2)=11
!      mt2ciw(5,3,2)=12


!      allocate(corr_iter(0:6,2))
!      allocate(at_bonding(0:6,2))
!      corr_iter=0
!      corr_iter(0,2)=2
!      corr_iter(1,1)=2
!      corr_iter(1,2)=2
!      corr_iter(2,2)=2
!      at_bonding=0


!      num_orb_cor=12
!      num_orb_iter=0
!      num_orb_proj=num_orb_cor+num_orb_iter
!      allocate(at_orb_ciw(2,num_orb_proj))
!      at_orb_ciw(1,1)=1
!      at_orb_ciw(2,1)=1
!      at_orb_ciw(1,2)=1
!      at_orb_ciw(2,2)=2
!      at_orb_ciw(1,3)=1
!      at_orb_ciw(2,3)=3
!      at_orb_ciw(1,4)=2
!      at_orb_ciw(2,4)=1
!      at_orb_ciw(1,5)=2
!      at_orb_ciw(2,5)=2
!      at_orb_ciw(1,6)=2
!      at_orb_ciw(2,6)=3
!      at_orb_ciw(1,7)=2
!      at_orb_ciw(2,7)=4
!      at_orb_ciw(1,8)=2
!      at_orb_ciw(2,8)=5
!      at_orb_ciw(1,9)=2
!      at_orb_ciw(2,9)=6
!      at_orb_ciw(1,10)=2
!      at_orb_ciw(2,10)=7
!      at_orb_ciw(1,11)=2
!      at_orb_ciw(2,11)=8
!      at_orb_ciw(1,12)=2
!      at_orb_ciw(2,12)=9
!      natom_ciw=2
!      n_ciw_max=9
!      allocate(n_ciw_atom(2))
!      n_ciw_atom(1)=3
!      n_ciw_atom(2)=9
!      allocate(iatom_ciw(2))
!      iatom_ciw(1)=1
!      iatom_ciw(2)=2

!      n_ciw_max=maxval(n_ciw_atom)
!      allocate(ciw_ind(n_ciw_max,natom_ciw))
!      ciw_ind(1,1)= 1
!      ciw_ind(2,1)=2
!      ciw_ind(3,1)=3
!      ciw_ind(1,2)=4
!      ciw_ind(2,2)=5
!      ciw_ind(3,2)=6
!      ciw_ind(4,2)=7
!      ciw_ind(5,2)=8
!      ciw_ind(6,2)=9
!      ciw_ind(7,2)=10
!      ciw_ind(8,2)=11
!      ciw_ind(9,2)=12


      r1_begin(1)=0.0
      r1_begin(2)=0.0
      r1_begin(3)=0.0
      r1_end(1)=0.0
      r1_end(2)=0.0
      r1_end(3)=1.0
      ndiv_r1 =100
!      enpj_emin=-20.0*2.d0/evolt
!      enpj_emax=20.0*2.d0/evolt
!      u_mode='wnds'
!      p_low_mode = 1
!      p_local_r = .false.
!      p_local_k = .false.
!      u_linear = .false.
      we_need_u = .true.
      w_q_zero = 0
      nbndf_bnd_0(1)=0
      nbndf_bnd_0(2)=0
!      projemax=80.0*2.d0/evolt
!      projemin=-80.0*2.d0/evolt
      n_pb_s_max = 80
      eps_pb_vrt=1.e-03
      eps_v_coul=1.e-04
      p_local_r =.false.
      v_repr = 'NON'
!      write(100,*)'eps_pb_vrt',eps_pb_vrt,
!     $ 'eps_v_coul',eps_v_coul 
!      ewald = 2.00
!      ewald0 = 2.00
      obietto='   solid'
      mem_using=0.d0
      mem_max=0.d0
      mem_where='qft'

!      r_ref=0.d0
!      nproc_tau=24
!      nproc_k=16
!      n_tau=46
!      n_nu_exa=29
!      n_nu_asy=18
!      n_nu=n_nu_exa+n_nu_asy
!      n_omega_exa=29
!      n_omega_asy=18
!      n_omega=n_omega_exa+n_omega_asy



      call read_ini
#ifdef USE_HDF5
      call read_ntau_uca_info_hdf5
#else
      call read_ntau_uca_info
#endif
      if(ubi=='dft') then
       ubi_num=1
      else if(ubi==' qp') then
       ubi_num=4
      endif      
      
c$$$      if (me .lt. 10) then
c$$$        write(ch, '(a, i1, a)') "check_", me, ".dat"
c$$$      elseif (me .lt. 100) then
c$$$        write(ch, '(a, i2, a)') "check_", me, ".dat"
c$$$      elseif (me .lt. 1000) then
c$$$        write(ch, '(a, i3, a)') "check_", me, ".dat"
c$$$      elseif (me .lt. 10000) then
c$$$        write(ch, '(a, i4, a)') "check_", me, ".dat"
c$$$      endif
c$$$
c      if (me .lt. 10) then
c        write(chbk, '(a, i1, a)') "check_", me, ".bkt"
c      elseif (me .lt. 100) then
c        write(chbk, '(a, i2, a)') "check_", me, ".bkt"
c      elseif (me .lt. 1000) then
c        write(chbk, '(a, i3, a)') "check_", me, ".bkt"
c      elseif (me .lt. 10000) then
c        write(chbk, '(a, i4, a)') "check_", me, ".bkt"
c      endif


c$$$      open(179, file=trim(ch))
c      write(179,*)'bkk-test',me
c      open(100, file=trim(chbk))
      iun=993
      if(maswrk) then
        open(iun,file='OUT')
!        write(iun,'(a)')  "git info"
!        write(iun,'(a)')  gitdate
!        write(iun,'(a)')  gitbranch
        write(iun,*)  
        write(iun,*)  
      endif


      if(maswrk) write(iun,*)' Number of processors is ',nproc

 
      call start_par_mb
c      if (me3_k .lt. 10) then
c        write(me3_k_ch, '(i1)') me3_k
c      elseif (me3_k .lt. 100) then
c        write(me3_k_ch, '(i2)') me3_k
c      elseif (me3_k .lt. 1000) then
c        write(me3_k_ch, '(i3)') me3_k
c      elseif (me3_k .lt. 10000) then
c        write(me3_k_ch, '(i4)') me3_k
c      endif
c      if (me3_tau .lt. 10) then
c        write(me3_tau_ch, '(i1)') me3_tau
c      elseif (me3_tau .lt. 100) then
c        write(me3_tau_ch, '(i2)') me3_tau
c      elseif (me3_tau .lt. 1000) then
c        write(me3_tau_ch, '(i3)') me3_tau
c      elseif (me3_tau .lt. 10000) then
c        write(me3_tau_ch, '(i4)') me3_tau
c      endif
c      if (me .lt. 10) then
c        write(me_ch, '(i1)') me
c      elseif (me .lt. 100) then
c        write(me_ch, '(i2)') me
c      elseif (me .lt. 1000) then
c        write(me_ch, '(i3)') me
c      elseif (me .lt. 10000) then
c        write(me_ch, '(i4)') me
c      else 
c        write(me_ch, '(i5)') me
c      endif
c      write(179, '(3(a, i5, 2x))')
c     $  'me3_k', me3_k, 'me3_tau', me3_tau, 'me', me

!      call start_par_solid
#ifdef USE_HDF5
      call read_uca_info_hdf5
#else
      call read_uca_info
#endif
      call get_rr
!      call start_par_mbi
      iun=993
      call read_ini
      if(ubi=='dft') then
        if(allocated(q_qp)) deallocate(q_qp)
        allocate(q_qp(nbndf,nbndf,ndim3_k(me3_k+1),nspin_0))
        e_qp=e_bnd
        z_qp=z_bnd
        ev_qp=ev_bnd
        q_qp=(0.d0,0.d0)
        chem_pot_qp=chem_pot
        do ispin=1,nspin
          do ind_k=1,ndim3_k(me3_k+1)
            k=n3_mpi_k(me3_k+1)+ind_k
            do i=1,n_bnd(k,ispin)
              q_qp(i,i,ind_k,ispin)=dcmplx(1.d0,0.d0) !assuming that calculation is from hartree fock
            enddo
          enddo
        enddo
      endif
      itest=1
      v_mm=(0.d0,0.d0)

      if(allocated(phase_k_rr)) then 
         deallocate(phase_k_ rr)
      endif
	  na=max(ndiv(1)-1,3)
	  nb=max(ndiv(2)-1,3)
	  nc=max(ndiv(3)-1,3)
        allocate(phase_k_rr(nqdiv,-na:na,-nb:nb,-nc:nc))	
	  pi2=pi+pi
	  do ic=-nc,nc
	    do ib=-nb,nb
	      do ia=-na,na
              v(:)=ia*rbas(:,1)+ib*rbas(:,2)+ic*rbas(:,3)
              do k=1,nqdiv
	          fas=pi2*dot_product(pnt(:,k),v)
	          phase_k_rr(k,ia,ib,ic)=dcmplx(cos(fas),sin(fas))
	        enddo
	      enddo
	    enddo
	  enddo


!       if (me3_tau .eq. 0) then
!        do ind_k=1,ndim3_k(me3_k+1)
!      do ispin=1, nspin_0
!      write(100,*)'q_qp',q_qp(:,:,ind_k,ispin)
!      enddo
!      enddo
!      endif

      ubi='dmf'
      ubi_num=7
      call mpi_bcast(nproc_k,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      call cal_weig_rc
!      seed_name='wannier'
      call cor_and_iter_wan
      call mpi_bcast(nproc_k,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
!      write(100,*)"ev_bnd",ev_bnd,"g_full",g_full,
!     $  "ab_nu_from_tau",ab_nu_from_tau,"v_ii",v_ii

!        write(100,*)"dr",dr,"eps_pb",eps_pb,"fun_pb",fun_pb,
!     $ "gfun",gfun,"gfund",gfund,"h",h,"ind_pb",ind_pb,
!     $ "ind_prod",ind_prod,"ind_prod0",ind_prod0,"ind_wf",ind_wf,
!     $ "iopb",iopb,"iopb0",iopb0,"irel",irel,"is",is,"iun",iun,
!     $ "lf_pb",lf_pb,"lfun",lfun_pb,"li_val",li_val,
!     $ "lim_pb_mt",lim_pb_mt,"lm_pbmt",lm_pbmt,"lmpb",lmpb,
!     $ "lval_pb",lval_pb,"max_e_pb",max_e_pb,"maxlfpb",maxlfpb,
!     $ "maxnrad",maxnrad,"maxpb",maxpb,"maxpbr",maxpbr,
!     $ "n_pbmt",n_pbmt,"n_pbmt0",n_pbmt0,"n_pbmtm",n_pbmtm,
!     $ "natom",natom,"nhsort",nhsort,"nrad",nrad,"nsort",nsort,
!     $ "nspin",nspin,"ntle_pb",ntle_pb,"r",r,"txtel",txtel

!        write(100,*)"e_qp",e_qp,"include_bands",include_bands,
!     $ "num_bands",num_bands,
!     $ "chem_pot",chem_pot,"evolt",evolt,"num_bands",num_bands

!      clight=274.074e+00
!       write(100,*)"z_bnd",z_bnd

! important
!       e_qp=e_bnd
!       z_qp=z_bnd
!       ev_qp=ev_bnd
!       q_qp=(0.d0,0.d0) 
!       chem_pot_qp=chem_pot
!       do ispin=1,nspin
!         do ind_k=1,ndim3_k(me3_k+1)
!           k=n3_mpi_k(me3_k+1)+ind_k
!           do i=1,n_bnd(k,ispin)
!             q_qp(i,i,ind_k,ispin)=dcmplx(1.d0,0.d0) !assuming that calculation is from hartree fock
!           enddo
!         enddo
!       enddo
      call solid_scf_new
c      close(179)

      if(maswrk) then
        close(iun)

      endif
      call ending_qft
      end
