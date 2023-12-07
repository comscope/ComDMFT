      Program comlowh

      use comlowh_mod

      implicit none
      include "mpif.h"
      integer ::
     $  ii,jj,kk,irk,iorb,jorb,iatom,icell,iio,iomega,tempk(3)
      double precision ::pi,vol_ratio,tempdble3(3)
      character*30 :: ch
      double precision, allocatable :: 
     $  tempdble(:,:,:,:)
      complex*16 :: ai

      CALL MPI_INIT(mpierr)
      
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,MPIERR)
      ME = irank
      me1=me+1
      MASTER = 0                                                        
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,MPIERR)      


      NPROC = isize
      MASWRK = ME .EQ. MASTER

      iun=993

      pi=datan2(1.0d0, 1.0d0)*4.0d0
      ai=dcmplx(0.0d0, 1.0d0)

      open(unit=iun, file='comlowh.log')

c$$$      if (me .lt. 10) then
c$$$        write(ch, '(a, i1, a)') "check_comlowh_", me, ".dat"
c$$$      elseif (me .lt. 100) then
c$$$        write(ch, '(a, i2, a)') "check_comlowh_", me, ".dat"
c$$$      elseif (me .lt. 1000) then
c$$$        write(ch, '(a, i3, a)') "check_comlowh_", me, ".dat"
c$$$      elseif (me .lt. 10000) then
c$$$        write(ch, '(a, i4, a)') "check_comlowh_", me, ".dat"
c$$$      elseif (me .lt. 100000) then
c$$$        write(ch, '(a, i5, a)') "check_comlowh_", me, ".dat"
c$$$      elseif (me .lt. 1000000) then
c$$$        write(ch, '(a, i6, a)') "check_comlowh_", me, ".dat"        
c$$$      endif

c$$$      open(179, file=trim(ch))

c$$$      write(179,*) 'here'
      seed_name='wannier'      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     1. read parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! the next block reads in comlowh.ini
      if (maswrk) then ! do this only on the root node

        open(unit=10, file='comlowh.ini', status='old')

        read(10,*) cal_mode
! from Sang:  dos_qp and band are before analytic continuation, and dos and spectral are after analytic continuation.
! band shows quasiparticle bands, and can plot the results by running: python ./band.py
! spectral shows A(k,omega), and can plot the results of spectral by running:  python band_spectra.py

! the only references to cal_mode=4,5 [no analytic continuation] are in the next if/else statement.
!     cal_mode 1 => delta
!       This not postprocessing.

!     'dos_qp' => cal_mode 4 => delta imag qp hamiltonian
!       This postprocessing shows the DOS without doing analytic continuation.

!     'band' => cal_mode 5 => delta imag path qp hamiltonian
!       This postprocessing shows quasiparticle bands without doing analytic continuation.

!     'dos' => cal_mode 2 => delta real
!       This postprocessing shows the DOS after analytic continuation.


!     'spectral' => cal_mode 3 => delta real path
!       This postprocessing shows A(k,omega) .

        if (cal_mode .eq. 1) then
          is_matsubara=1 ! frequency(iomega) is imaginary
          is_kpath=0
          is_recal_ef=1 !  if is_metal=0 then is_recal_ef will be set to zero.  if is_recal_ef=1 then we will recalculate the Fermi energy and save it in ef.dat.  Otherwise we will read in the Fermi energy from ef.dat.
! Is ef.dat is set to zero by comdmft.py/prepare_initial_ef before the first call to ComLowH in self-consistent dft+dmft or qsgw+dmft.
! If doing self-consistent dft+dmft or qsgw+dmft, AND if not(is_metal==0), ComLowH never reads in ef.dat; it only recomputes it and saves it out during each iteration.
! If (is_metal==0), or if doing post-processing, ComLowH reads in ef.dat and uses its value.
! ef.dat is not used by other parts of ComDMFT.
          is_qp=0 ! use cal_g not cal_h_qp
        elseif (cal_mode .eq. 2) then ! yes needs analytic continuation, DOS
          is_matsubara=0 ! frequency(iomega) is real
          is_kpath=0
          is_recal_ef=0
          is_qp=0 ! use cal_g not cal_h_qp
        elseif (cal_mode .eq. 3) then ! yes needs analytic continuation, A(k,omega)
          is_matsubara=0 ! frequency(iomega) is real
          is_kpath=1 ! is_kpath means we are calculating either A((k,omega) [if cal_mode=3] or qp spectra [if cal_mode=5]
          is_recal_ef=0
          is_qp=0 ! use cal_g not cal_h_qp
        elseif (cal_mode .eq. 4) then ! the output files are exactly the same as (cal_mode .eq. 5)
          is_matsubara=1 ! frequency(iomega) is imaginary
          is_kpath=0
          is_recal_ef=0
          is_qp=1 ! use cal_h_qp
        elseif (cal_mode .eq. 5) then ! no analytic continuation, produces qp spectra, does not produce a DOS or A(k,omega)
          is_matsubara=1 ! frequency(iomega) is imaginary
          is_kpath=1 ! is_kpath means we are calculating either A(k,omega) [if cal_mode=3] or qp spectra [if cal_mode=5]
          is_recal_ef=0
          is_qp=1 ! ! use cal_h_qp
        endif
        
        num_cell=1

c$$$        read(10,*) num_cell
        read(10,*) num_atom ! len(control['impurity_wan'])
        allocate(num_orb(num_atom))

        num_orb=0
        read(10,*) (num_orb(iatom), iatom=1, num_atom) ! nimp_orb
        num_orb_max=maxval(num_orb)

        allocate(iac(num_atom,num_cell))

# read in map(str,control['impurity_problem_equivalence'])
        read(10,*)
     $    ((iac(iatom,icell),iatom=1,num_atom),icell=1,num_cell)
        
        if (minval(iac) .lt. 0) then
          n_iac=maxval(iac)*2
          n_iac_nm=maxval(iac)
          n_iac_mat=n_iac+1
          n_iac_mat_i=-n_iac_nm
          n_iac_mat_f=n_iac_nm          
          is_magnetic=1
        else
          n_iac=maxval(iac)
          n_iac_nm=maxval(iac)
          n_iac_mat=n_iac
          n_iac_mat_i=1
          n_iac_mat_f=n_iac_nm                    
          is_magnetic=0
        endif

        allocate(ndeg_iac(n_iac_mat_i:n_iac_mat_f))
        ndeg_iac=0
        allocate(norb_iac(n_iac_mat_i:n_iac_mat_f))
        norb_iac=0
        
        do icell=1, num_cell ! num_cell=1
          do iatom=1, num_atom
            norb_iac(iac(iatom,icell))=num_orb(iatom)
            ndeg_iac(iac(iatom,icell))
     $        =ndeg_iac(iac(iatom,icell))+1
          enddo
        enddo
        
        allocate(ioac(num_orb_max, num_orb_max,
     $    n_iac_mat_i:n_iac_mat_f))
        ioac=0          

# This reads  a matrix, from the entries of  control['impurity_problem_equivalence'].
        do ii=n_iac_mat_i, n_iac_mat_f
          if ((is_magnetic .eq. 1) .and. (ii .eq. 0)) cycle
          do iorb=1, norb_iac(ii)
            read(10,*) (ioac(iorb,jorb,ii),
     $        jorb=1, norb_iac(ii))
          enddo
        enddo
        n_ioac=maxval(ioac)        

        allocate(iiiio(n_ioac,n_iac_mat_i:n_iac_mat_f))
        allocate(iio_diagonal(n_ioac,n_iac_mat_i:n_iac_mat_f))        
        iiiio=0
        n_i=0
        do ii=n_iac_mat_i, n_iac_mat_f
          if ((is_magnetic .eq. 1) .and. (ii .eq. 0)) cycle          
          do iio=1, maxval(ioac(:,:,ii))
            n_i=n_i+1
            iiiio(iio,ii)=n_i
          enddo
          do iorb=1, norb_iac(ii)
            do jorb=1, norb_iac(ii)
              iio=ioac(iorb,jorb,ii)
              if (iio .eq. 0) cycle
              if (iorb .eq. jorb) then
                iio_diagonal(iio,ii)=1
              else
                iio_diagonal(iio,ii)=0
              endif
            enddo
          enddo
        enddo

        
        allocate(ndeg_ioac(n_ioac, n_iac_mat_i:n_iac_mat_f))
        ndeg_ioac=0
        allocate(ndeg_itot(n_ioac, n_iac_mat_i:n_iac_mat_f))
        ndeg_itot=0          
        
        
        do ii=n_iac_mat_i, n_iac_mat_f
          if ((is_magnetic .eq. 1) .and. (ii .eq. 0)) cycle
          do iorb=1, norb_iac(ii)
            do jorb=1, norb_iac(ii)
              iio=ioac(iorb,jorb,ii)
              if (iio .eq. 0) cycle
              ndeg_ioac(iio,ii)
     $          =ndeg_ioac(iio,ii)+1
            enddo
          enddo
        enddo
        
        do icell=1, num_cell ! num_cell=1
          do iatom=1, num_atom
            ii=iac(iatom,icell)              
            do iorb=1, num_orb(iatom)
              do jorb=1, num_orb(iatom)
                iio=ioac(iorb,jorb,ii)
                if (iio .eq. 0) cycle                
                ndeg_itot(iio,ii)
     $            =ndeg_itot(iio,ii)+1
              enddo
            enddo
          enddo
        enddo
        
        ndeg_ioac_max=maxval(ndeg_ioac)
        
        
        allocate(wan_correlated(num_orb_max, num_atom))
        wan_correlated=0

# Read in indices of the orbitals used for impurities, which came from control['impurity_wan'][iatom].
        read(10,*) ((wan_correlated(iorb,iatom),
     $    iorb=1,num_orb(iatom)),iatom=1,num_atom)

        read(10,*) ewin_min, ewin_max ! control['proj_win_min'], control['proj_win_max']
        read(10,*) n_omega ! control['n_omega']
        read(10,*) gamma_selfenergy ! 0.0
        read(10,*) gamma_real ! 0.0
        dc_shift=0.0d0
        read(10,*) beta ! imp['beta']
        read(10,*) delta_charge !! for one cell. system charge ! control['doping']
        read(10,*) is_metal ! from is_recal_ef in python
c$$$        is_metal=1
        is_recal_ef=is_recal_ef*is_metal
        read(10,*) den_mode ! 'bnd'
        read(10,*) is_spinorbit ! control['spin_orbit']

! ndiv and nqdiv are the coarse grid, ndim_kk_co,n_mpi_kk_co
! ndiv_fi and nqdiv_fi are the fine grid, ndim_kk_fi,n_mpi_kk_fi

        !  if (is_kpath .eq. 0), comlowh uses the ndiv_fi to set up kpt_latt_fi, and then writes it out to the kpoints file.
        ! if not(is_kpath .eq. 0), then comlowh reads from the kpoints file into kpt_latt_fi
        ! is_kpath means we are calculating either A((k,omega) [if cal_mode=3] or qp spectra [if cal_mode=5]
        ! the comdmft.ini file can contain a kpoints section, which may be a list of kpoints or may instead be reference to a kpoints file.
        if (is_kpath .eq. 0) then
          read(10,*) ndiv_fi ! wan_hmat['kgrid']
          nqdiv_fi=ndiv_fi(1)*ndiv_fi(2)*ndiv_fi(3)
          allocate(kpt_latt_fi(3, nqdiv_fi))
          kpt_latt_fi=0.0d0
          
          do irk=1, nqdiv_fi
            call indexing(nqdiv_fi,3,ndiv_fi,0,irk,tempk)
            do ii=1, 3
              kpt_latt_fi(ii,irk)
     $          =dble(tempk(ii)-1.0d0)/dble(ndiv_fi(ii))
            enddo
          enddo
          open(unit=11, file='kpoints')
          do irk=1, nqdiv_fi          
            write(11, '(i5,2x,3(f20.12,2x))') irk, kpt_latt_fi(:,irk)
          enddo
          close(11)
        else
          open(unit=11, file='kpoints', status='old')
          read(11,*)
          read(11,*)           
          read(11,*) nqdiv_fi
          if (nqdiv_fi .gt. 0) then
            allocate(kpt_latt_fi(3,nqdiv_fi))
            kpt_latt_fi=0.0d0
            do ii=1, nqdiv_fi
              read(11,*) jj, (kpt_latt_fi(jj, ii), jj=1,3)
            enddo
          endif
          close(11)
        endif ! Done setting up the kpoints file.

        close(10)
! this ends the block that reads in comlowh.ini


c$$$  if (num_cell .ne. 1) then !! in angstorm
c$$$  allocate(lpoint(3,num_cell))
c$$$  lpoint=0.0d0
c$$$  read(10,*) (real_lattice_conv(ii, 1), ii=1,3)
c$$$  read(10,*) (real_lattice_conv(ii, 2), ii=1,3)
c$$$  read(10,*) (real_lattice_conv(ii, 3), ii=1,3)
c$$$  do icell=1, num_cell
c$$$  read(10,*) (lpoint(ii, icell), ii=1,3)
c$$$  enddo
c$$$  endif
#ifdef USE_HDF5
#else

        ! this is a file read
        ! reads in wannier.dat
        call read_wann_dat
#endif

      endif ! done reading in ini files, only on the root node
      




!            # comlowh also writes to other files, some of which are for postprocessing
!            # some of the files are:
!            #   cal_mode=2: tdos.dat, pdos.dat, inquires about momentum_optics_wan.dat
!            #   cal_mode=3: spectral.matrix, spectral.dat, spectral_orb.dat, wannier_band_non_interpolated.dat

      ! This next block broadcasts ini file info from the root node to other nodes
      call mpi_bcast(cal_mode,1,mpi_integer,0,MPI_COMM_WORLD,mpierr) !
      call mpi_bcast(is_matsubara,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(is_kpath,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(is_recal_ef,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(is_qp,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)                        
      call mpi_bcast(num_cell,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)      
      call mpi_bcast(num_atom,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      if (.not. maswrk) allocate(num_orb(num_atom))      
      call mpi_bcast(num_orb,num_atom,mpi_integer,0,MPI_COMM_WORLD,
     $  mpierr)
      call mpi_bcast(num_orb_max,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      if(.not. maswrk) allocate(iac(num_atom,num_cell))
      call mpi_bcast(iac,num_atom*num_cell,mpi_integer,0,MPI_COMM_WORLD,
     $  mpierr)      
      call mpi_bcast(n_iac,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(n_iac_nm,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(n_iac_mat,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(n_iac_mat_i,1,mpi_integer,0,
     $  MPI_COMM_WORLD,mpierr)
      call mpi_bcast(n_iac_mat_f,1,mpi_integer,0,
     $  MPI_COMM_WORLD,mpierr)            
      call mpi_bcast(is_magnetic,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      if (.not. maswrk) allocate(norb_iac(n_iac_mat_i:n_iac_mat_f))
      call mpi_bcast(norb_iac,n_iac_mat,mpi_integer,
     $  0,mpi_comm_world,mpierr)
      if (.not. maswrk) allocate(ndeg_iac(n_iac_mat_i:n_iac_mat_f))
      call mpi_bcast(ndeg_iac,n_iac_mat,mpi_integer,
     $  0,mpi_comm_world,mpierr)
      if (.not. maswrk) allocate(ioac(num_orb_max, num_orb_max,
     $  n_iac_mat_i:n_iac_mat_f))
      call mpi_bcast(ioac,num_orb_max**2*n_iac_mat,mpi_integer,
     $  0,mpi_comm_world,mpierr)
      call mpi_bcast(n_ioac,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)      
      if (.not. maswrk) allocate(iiiio(n_ioac,n_iac_mat_i:n_iac_mat_f))
      call mpi_bcast(iiiio,n_iac_mat*n_ioac,mpi_integer,
     $  0,mpi_comm_world,mpierr)
      if (.not. maswrk)
     $  allocate(iio_diagonal(n_ioac,n_iac_mat_i:n_iac_mat_f))
      call mpi_bcast(iio_diagonal,n_iac_mat*n_ioac,mpi_integer,
     $  0,mpi_comm_world,mpierr)
      call mpi_bcast(n_i,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)       

      if (.not. maswrk)
     $  allocate(ndeg_ioac(n_ioac, n_iac_mat_i:n_iac_mat_f))              
      call mpi_bcast(ndeg_ioac,n_ioac*n_iac_mat,
     $  mpi_integer,0,MPI_COMM_WORLD,mpierr)
      if (.not. maswrk)      
     $  allocate(ndeg_itot(n_ioac, n_iac_mat_i:n_iac_mat_f))      
      call mpi_bcast(ndeg_itot,n_ioac*n_iac_mat,
     $  mpi_integer,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(ndeg_ioac_max,1,
     $  mpi_integer,0,MPI_COMM_WORLD,mpierr)
      

      if (.not. maswrk) allocate(wan_correlated(num_orb_max,num_atom))
      call mpi_bcast(wan_correlated,num_atom*num_orb_max,mpi_integer,
     $  0,MPI_COMM_WORLD,mpierr)
      
      call mpi_bcast(ewin_min,1,mpi_double_precision,
     $  0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(ewin_max,1,mpi_double_precision,
     $  0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(n_omega,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(gamma_selfenergy,1,mpi_double_precision,
     $  0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(gamma_real,1,mpi_double_precision,
     $  0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(dc_shift,1,mpi_double_precision,
     $  0,MPI_COMM_WORLD,mpierr)      
      
      call mpi_bcast(beta,1,mpi_double_precision,
     $  0,MPI_COMM_WORLD,mpierr)      
      call mpi_bcast(delta_charge,1,
     $  mpi_double_precision,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(is_metal,1,
     $  mpi_integer,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(den_mode,3,
     $  mpi_character,0,MPI_COMM_WORLD,mpierr)            
      call mpi_bcast(is_spinorbit,1,
     $  mpi_integer,0,MPI_COMM_WORLD,mpierr)
      ! is_kpath means we are calculating either A((k,omega) [if cal_mode=3] or qp spectra [if cal_mode=5]
      if (is_kpath .eq. 0) then      
        call mpi_bcast(ndiv_fi,3,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      endif
      call mpi_bcast(nqdiv_fi,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      if (.not. maswrk) allocate(kpt_latt_fi(3, nqdiv_fi))
      call mpi_bcast(kpt_latt_fi,3*nqdiv_fi,mpi_double_precision,0,
     $  MPI_COMM_WORLD,mpierr)      
#ifdef USE_HDF5
       call read_wannier_dat_hdf5
#else
      
      call mpi_bcast(real_lattice,9,mpi_double_precision,
     $  0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(recip_lattice,9,mpi_double_precision,
     $  0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(num_bands,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(num_wann,1,mpi_integer,
     $  0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(ndiv,3,mpi_integer,
     $  0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(nqdiv,1,mpi_integer,
     $  0,MPI_COMM_WORLD,mpierr)
      if (.not. maswrk) allocate(kpt_latt(3,nqdiv))
      call mpi_bcast(kpt_latt,3*nqdiv,mpi_double_precision,
     $  0,MPI_COMM_WORLD,mpierr)
      if (.not. maswrk) allocate(include_bands(num_bands))
      call mpi_bcast(include_bands,num_bands,mpi_integer,0,
     $  MPI_COMM_WORLD,mpierr)
      if (.not. maswrk) allocate(v_matrix(num_bands,num_wann,nqdiv))
      call mpi_bcast(v_matrix,num_bands*num_wann*nqdiv,
     $  mpi_double_complex,0,MPI_COMM_WORLD,mpierr)
      
      if (.not. maswrk) allocate(eigenvalues(num_bands,nqdiv))
      call mpi_bcast(eigenvalues,num_bands*nqdiv,
     $  mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

      call mpi_bcast(nrpts,1,mpi_integer,0,
     $  MPI_COMM_WORLD,mpierr)
      if (.not. maswrk) allocate(irvec(3,nrpts))
      call mpi_bcast(irvec,3*nrpts,mpi_integer,0,
     $  MPI_COMM_WORLD,mpierr)      
      
      if (.not. maswrk) allocate(ndegen(nrpts))
      call mpi_bcast(ndegen,nrpts,mpi_integer,0,
     $  MPI_COMM_WORLD,mpierr)

      call mpi_bcast(nelectron,1,
     $  mpi_double_precision,0,MPI_COMM_WORLD,mpierr)
#endif
! End of the block that broadcasts ini file info from the root node to other nodes.


      vol_ratio=1.0d0

      
c$$$      if (num_cell .ne. 1) then
c$$$        call mpi_bcast(real_lattice_conv,9,
c$$$     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)
c$$$        if (.not. maswrk) allocate(lpoint(3,num_cell))
c$$$        call mpi_bcast(lpoint,3*num_cell,
c$$$     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)
c$$$      endif

      num_proj=sum(num_orb)*num_cell !  ! num_cell=1
      num_wancell=num_wann*num_cell ! num_wancell = num_wann
      num_rest=num_bands-num_wann

      allocate(wan_cor_index(num_orb_max,num_atom,num_cell))

! wan_cor_index = wan_correlated
! for instance, for a single d-shell atom, wan_cor_index=wan_correlated might be [8,9,10,11,12], and its iatom and icell indices are constrained to 1.
      wan_cor_index=0
      do icell=1, num_cell !  ! num_cell=1
        do iatom=1, num_atom
          do iorb=1, num_orb(iatom)
            wan_cor_index(iorb,iatom,icell)
     $        =wan_correlated(iorb,iatom)+(icell-1)*num_wann
          enddo
        enddo
      enddo

c$$$      write(179,*) 'here'
      allocate(frequency_dble(n_omega)) !
      
      ef=0.0d0
      allocate(sig(n_ioac,n_iac_mat_i:n_iac_mat_f,n_omega))
      allocate(dc(n_ioac,n_iac_mat_i:n_iac_mat_f))        
      sig=0.0d0
      dc=0.0d0
      


! This reads sig.dat into sig.
! Here is is the origin of sig.dat:
! comdmft.py's run_impurity_solver runs CTQMC and EVALSIM, and then:
! reads in ctqmc's output from params.obs.json and updates convergence.log
! "self-energy" from params.obs.json is saved in sig_bare.dat
! "self-energy" is smoothed using Gaussian broadening and stored in sigma.
!    sigma is saved in  sig_smth.dat.
! If any element of imag(sigma) is positive, then sig_causality is set
!       to 0=False; otherwise it is 1=True.
! If any element of imag(sigma) is positive, then sigma_to_delta = sigma_old [read in from sig.dat].
!    If it is true then sigma_to_delta is a mix of sigma with sigma_old [ read in from sig.dat].
! sigma_to_delta is saved in sig.dat.
! In summary, sig.dat comes from params.obs.json, then gaussian smoothing, then a check for negative imaginary parts, then mixing.
! ------------
! On the initial self-consistent cycle of dft+dmft and lqsgw+dmft, things are different:
! comdmft.py's generate_initial_self_energy creates sig.dat.
! If ('initial_self_energy' in control,  comdmft.py copies that self_energy to sig.dat.
!   Possibly it also copies information from initial_impurity_dir.
! Otherwise, comdmft.py copies data from dc.dat to sig.dat.
! ----------------
! If doing lqsgw+dmft, then before doing the self-consistent cycle there is a comlowh+comdc run.  To prepare for this dc.dat is all zeros, and sig.dat contains all zero's and omega's.  This is set up by comdmft.py's prepare_seed_dc_sig_and_wannier_dat.
      if (maswrk) then
        open(unit=10, file='sig.dat', status='old')
        allocate(tempdble(2,n_i,1,1))
        read(10,*) 
        do iomega=1, n_omega
          tempdble=0.0d0
          read(10,*) frequency_dble(iomega),
     $      ((tempdble(jj,ii,1,1), jj=1, 2), 
     $      ii=1, n_i)
          do ii=n_iac_mat_i, n_iac_mat_f
            if ((is_magnetic .eq. 1) .and. (ii .eq. 0)) cycle
            do iio=1, maxval(ioac(:,:,ii))
              sig(iio,ii,iomega)=
     $          dcmplx(tempdble(1,iiiio(iio,ii),1,1),
     $          tempdble(2,iiiio(iio,ii),1,1))
            enddo
          enddo
        enddo
        deallocate(tempdble)
        close(10)
      endif
      
      call mpi_bcast(sig,n_ioac*n_iac_mat*n_omega,
     $  mpi_double_complex,0,MPI_COMM_WORLD,mpierr)
      
      call mpi_bcast(frequency_dble,n_omega,
     $  mpi_double_precision,0,MPI_COMM_WORLD,mpierr)
      temperature=frequency_dble(1)/pi/(8.6173303d-5)
      if (maswrk) write(iun, *) 'temperature', temperature
      
! This reads dc.dat into dc .
! This is the origin of dc.dat: If doing qsgw+dmft, then run_dc runs ComDC, and then takes data (dc[0,jj,kk]) from ComDC's sig_mat.dat and puts it dc_mat.dat.  Alternatively, if doing dft+dmft, then run_dc/cal_nominal_dc creates dc_mat.dat.  Alternatively, if ('dc_mat_to_read' in control), then run_dc only copies an old version of dc_mat.dat from somewhere and uses that. After ComLowH is run, delta_postprocessing/cal_dc_diagonal takes takes from dc_mat.dat and puts it in dc.dat.
! when doing lqsgw+dmft, during the first call to comlowh + comdc (before the loop starts), dc.dat is all zeros.  This is set up by comdmft.py's prepare_seed_dc_sig_and_wannier_dat.
! dc.dat is used both by ComLowH and by CTQMC.
      if (maswrk) then
        open(unit=10, file='dc.dat', status='old')
        allocate(tempdble(2,n_i,1,1))
        tempdble=0.0d0
        read(10,*) ((tempdble(jj,ii,1,1), jj=1, 2), 
     $    ii=1, n_i)
        do ii=n_iac_mat_i, n_iac_mat_f
          if ((is_magnetic .eq. 1) .and. (ii .eq. 0)) cycle
          do iio=1, maxval(ioac(:,:,ii))
            dc(iio,ii)=
     $        dcmplx(tempdble(1,iiiio(iio,ii),1,1),
     $        tempdble(2,iiiio(iio,ii),1,1))
          enddo
        enddo
        deallocate(tempdble)
        close(10)

! This reads zinv_m1.dat into zinv_m1 .
! This is the origin of zinv_m1.dat: run_dc runs ComDC, and then takes data from ComDC's sig_mat.dat and puts it zinv_m1_mat.  After ComLowH is run, delta_postprocessing/cal_zinv_m1_diagonal takes takes info from zinv_m1_mat.dat and puts it in zinv_m1.dat
        inquire(file="zinv_m1.dat", exist=zinv_m1_exist)
        if (zinv_m1_exist)
     $    allocate(zinv_m1(n_ioac,n_iac_mat_i:n_iac_mat_f))        
        if (zinv_m1_exist) then
          open(unit=10, file='zinv_m1.dat', status='old')
          allocate(tempdble(2,n_i,1,1))
          tempdble=0.0d0
          read(10,*) ((tempdble(jj,ii,1,1), jj=1, 2), 
     $      ii=1, n_i)
          do ii=n_iac_mat_i, n_iac_mat_f
            if ((is_magnetic .eq. 1) .and. (ii .eq. 0)) cycle
            do iio=1, maxval(ioac(:,:,ii))
              zinv_m1(iio,ii)=
     $          dcmplx(tempdble(1,iiiio(iio,ii),1,1),
     $          tempdble(2,iiiio(iio,ii),1,1))
            enddo
          enddo
          deallocate(tempdble)
          close(10)
        endif
      endif

      call mpi_bcast(dc,n_ioac*n_iac_mat,
     $  mpi_double_complex,0,MPI_COMM_WORLD,mpierr)
      
      call mpi_bcast(zinv_m1_exist,1,
     $  mpi_logical,0,MPI_COMM_WORLD,mpierr)

      if (zinv_m1_exist) then
        if (.not. maswrk) then
          allocate(zinv_m1(n_ioac,n_iac_mat_i:n_iac_mat_f))
        endif
        call mpi_bcast(zinv_m1,n_ioac*n_iac_mat,
     $    mpi_double_complex,0,MPI_COMM_WORLD,mpierr)
      endif
      

      
      
c$$$      write(179,*) 'here'
      
      allocate(frequency(n_omega))
      frequency=0.0d0
c$$$      write(179,'(a, 1000(2f12.6,2x))') 'sig', sig(1,1,:)
c$$$      write(179,'(a, 1000(2f12.6,2x))') 'dc', dc(1,1)      
      

      do iomega=1,n_omega
        if (is_matsubara .eq. 1) then
          frequency(iomega)=dcmplx(0.0d0, frequency_dble(iomega))
        else
          frequency(iomega)=dcmplx(frequency_dble(iomega), gamma_real)
        endif
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     3. intermediate cell setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c$$$      if (num_cell .ne. 1) then
c$$$        call reciprocal_lattice_cal(real_lattice_conv,
c$$$     $    recip_lattice_conv)
c$$$        call vol_ratio_cal(real_lattice,
c$$$     $    real_lattice_conv, vol_ratio)
c$$$        write(179,*) 'vol_ratio', vol_ratio
c$$$        if (nint(vol_ratio) .ne. num_cell) then
c$$$          if (maswrk) then
c$$$            write(iun,*) "real_lattice_conv is wrong"
c$$$          endif
c$$$          CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
c$$$          CALL MPI_FINALIZE(mpierr)        
c$$$        endif
c$$$      endif
      
c$$$      write(179,*) 'here'
c$$$      write(179,*) 'from comlowh.ini'
c$$$      
c$$$      write(179,*) "correlated orbital information"
c$$$      write(179,'(a, 100i5)') 'cal_mode', cal_mode
c$$$      write(179,'(a, 100i5)') 'is_matsubara', is_matsubara
c$$$      write(179,'(a, 100i5)') 'is_kpath', is_kpath
c$$$      write(179,'(a, 100i5)') 'is_recal_ef', is_recal_ef
c$$$      write(179,'(a, 100i5)') 'num_cell', num_cell      
c$$$      write(179,'(a, 100i5)') 'num_atom', num_atom
c$$$      write(179,'(a, 100i5)') 'num_orb', (num_orb(iorb),iorb=1,num_atom)
c$$$      write(179,'(a, 100i5)') 'iac', iac
c$$$      write(179,'(a, 100i5)') 'n_iac', n_iac
c$$$      write(179,'(a, 100i5)') 'n_iac_nm', n_iac_nm
c$$$      write(179,'(a, 100i5)') 'n_iac_mat', n_iac_mat
c$$$      write(179,'(a, 100i5)') 'n_iac_mat_i', n_iac_mat_i
c$$$      write(179,'(a, 100i5)') 'n_iac_mat_f', n_iac_mat_f
c$$$      write(179,'(a, 100i5)') 'is_magnetic', is_magnetic
c$$$      write(179,'(a, 100i5)') 'n_ioac', n_ioac
c$$$      write(179,'(a, 100i5)') 'ndeg_iac', ndeg_iac
c$$$      write(179,'(a, 100i5)') 'norb_iac', norb_iac
c$$$      write(179,'(a, 1000i5)') 'ioac', ioac
c$$$      write(179,'(a, 1000i5)') 'iiiio', iiiio
c$$$      write(179,'(a, 1000i5)') 'ndeg_ioac_max', ndeg_ioac_max            
c$$$      write(179,'(a, 1000i5)') 'ndeg_ioac', ndeg_ioac      
c$$$      write(179,'(a, 1000i5)') 'ndeg_itot', ndeg_itot
c$$$      write(179,*) num_proj
c$$$      write(179,*) num_wann
c$$$      write(179,*) num_wancell
c$$$      write(179,'(100i5)')
c$$$     $  ((wan_correlated(ii,jj),ii=1,num_orb(jj)),jj=1,num_atom)
c$$$      write(179,'(100i5)')
c$$$     $  (((wan_cor_index(ii,jj,kk),ii=1,num_orb(jj)),
c$$$     $  jj=1,num_atom),kk=1,num_cell)      
c$$$      write(179,*) ewin_min, ewin_max
c$$$      write(179,*) "primitive cell information"
c$$$      write(179,'(9f12.6)') ((real_lattice(ii,jj),ii=1,3),jj=1,3) ! Real lattice
c$$$      write(179,'(9f12.6)') ((recip_lattice(ii,jj),ii=1,3),jj=1,3) ! Reciprocal lattice
c$$$      write(179,*) num_bands    ! Number of bands
c$$$      write(179,*) nqdiv        ! K-points
c$$$      if (cal_mode .eq. 1) then
c$$$        write(179,'(3i5)') (ndiv(ii),ii=1,3) ! M-P grid
c$$$      endif
c$$$      write(179,*) num_wann     ! num_wann
c$$$      write(179,*) nelectron*num_cell
c$$$c$$$      if (num_cell .gt. 1) then
c$$$c$$$        write(179,*) "conventional cell w.r.t. primitive cell vector"
c$$$c$$$        write(179,'(9f12.6)') ((real_lattice_conv(i,j),i=1,3),j=1,3) ! Real lattice
c$$$c$$$        write(179,'(9f12.6)') ((recip_lattice_conv(i,j),i=1,3),j=1,3) ! Reciprocal lattice
c$$$c$$$        write(179,*) "lpoint"
c$$$c$$$        do icell=1, num_cell
c$$$c$$$          tempdble3=matmul(transpose(recip_lattice)/2.0/pi,
c$$$c$$$     $      lpoint(:,icell))
c$$$c$$$          write(179,'(9f12.6)') (lpoint(i,icell),i=1,3),
c$$$c$$$     $      (tempdble3(i),i=1,3) ! Real lattice
c$$$c$$$        enddo
c$$$c$$$      endif
c$$$      write(179,*) "fine_cell information"
c$$$      write(179,*) 'ndiv_fi', ndiv_fi
c$$$      write(179,*) 'nqdiv_fi', nqdiv_fi
c$$$      write(179,*) "other details"
c$$$      write(179,*) 'n_omega', n_omega
c$$$      write(179,*) 'gamma_selfenergy', gamma_selfenergy
c$$$      write(179,*) 'gamma_real', gamma_real
c$$$      write(179,*) 'beta', beta
c$$$      write(179,*) 'delta_charge', delta_charge
c$$$      write(179,*) 'is_metal', is_metal
c$$$      write(179,*) 'is_spinorbit', is_spinorbit
c$$$      write(179,*) 'include_bands', include_bands
c$$$
c$$$      write(179,*) 'eigenvalues', eigenvalues(:,1)
c$$$      write(179,*) 'v_matrix', sum(abs(v_matrix)**2, 1)

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     5. parallelization setup
!!!!!!!!!!!!!!!!!!!!!!p!!!!!!!!!!!!!!!!!!!!!!

! ndim_kk_co and n_mpi_kk_co are a parallelization scheme:
!   It divides the nqdiv k-points across nproc processors.
      if (cal_mode .eq. 1) then ! doing a non-postprocessing run
        allocate(ndim_kk_co(nproc)) 
        ndim_kk_co=0            !
        allocate(n_mpi_kk_co(nproc)) !
        n_mpi_kk_co=0      
        call size_shift_par(nqdiv,nproc,ndim_kk_co,n_mpi_kk_co)
        ndim_kk_co_max=maxval(ndim_kk_co)            
      endif

        
! ndim_kk_fi and n_mpi_kk_fi are a parallelization scheme:
!   It divides the nqdiv_fi k-points across nproc processors.
      allocate(ndim_kk_fi(nproc)) 
      ndim_kk_fi=0              !
      allocate(n_mpi_kk_fi(nproc)) !
      n_mpi_kk_fi=0
      call size_shift_par(nqdiv_fi,nproc,ndim_kk_fi,n_mpi_kk_fi)
      ndim_kk_fi_max=maxval(ndim_kk_fi)


! This does comlowh's actual work.
      if (is_qp .eq. 0) then
        if (maswrk) then
          write(iun, *) 'cal_g run'
        endif
c$$$        write(179, *) 'cal_g run'        
        call cal_g(nqdiv_fi,ndim_kk_fi,n_mpi_kk_fi,ndim_kk_fi_max,
     $    kpt_latt_fi)
      else ! not (is_qp .eq. 0)
        if (maswrk) then
          write(iun, *) 'cal_h_qp run'
        endif 
c$$$        write(179, *) 'cal_h_qp run'        
        call cal_h_qp(nqdiv_fi,ndim_kk_fi,n_mpi_kk_fi,ndim_kk_fi_max,
     $    kpt_latt_fi)
      endif ! not (is_qp .eq. 0)


! deallocate 
      deallocate(num_orb)
      deallocate(iac)      
      deallocate(ndeg_iac)
      deallocate(norb_iac)
      deallocate(ioac)
      deallocate(iiiio)
      deallocate(iio_diagonal)            
      deallocate(ndeg_ioac)
      deallocate(ndeg_itot)
      deallocate(wan_correlated)
      deallocate(kpt_latt_fi)
      deallocate(kpt_latt)      
      deallocate(include_bands)
      deallocate(v_matrix)
      deallocate(eigenvalues)
      deallocate(irvec)
      deallocate(ndegen)      
      deallocate(wan_cor_index)
      deallocate(frequency_dble)
      deallocate(sig)
      deallocate(dc)
      deallocate(frequency)      
      if (cal_mode .eq. 1) then ! doing a non-postprocessing run
        deallocate(ndim_kk_co)
        deallocate(n_mpi_kk_co)
      endif
      deallocate(ndim_kk_fi)
      deallocate(n_mpi_kk_fi)

c$$$      if (num_cell .ne. 1) then !! in angstorm
c$$$        deallocate(lpoint)
c$$$      endif
c$$$      if (is_kpath .eq. 1) then
c$$$        deallocate(kpath)
c$$$      endif      
c$$$      close(179)


      if (maswrk) then
        close(iun)
      endif

      call mpi_finalize(mpierr)

      end  ! Program comlowh

      
      subroutine cal_g(numk,ndim_k,n_mpi_k,ndim_k_max,kvec)
! numk<=>nqdiv_fi,ndim_k<=>ndim_kk_fi,n_mpi_k<=>n_mpi_kk_fi, ndim_k_max<=>ndim_kk_fi_max,kvec<=>kpt_latt_fi
! ndim_k and n_mpi_k are a parallelization scheme for numk k-points distributed over nproc processors.  kvec is a complete list of all numk k-points. ndim_k_max=maxval(ndim_k)

      use comlowh_mod

      implicit none
      include "mpif.h"

      integer, intent(in) :: numk, ndim_k(nproc), n_mpi_k(nproc),
     $  ndim_k_max
      double precision, intent(in) :: kvec(3,numk)

      integer ::
     $  ii,jj,i1,j1,ik,irk,iomega,ibnd,jbnd,iwancell,jwancell,kwancell,
     $  iorb,jorb,iatom,icell,jcell,iwan,jwan,igshift,
     $  gvec_conv_shift(3,num_cell,numk),iio,fail_ind,proc_num,indout

      double precision ::pi,phaseval,phasetau,inv_numk,
     $  deltaij,diag(num_wancell)
      
      complex*16 :: ai, phaseval_dcmplx(num_cell),
     $  h_k_wan_co(nqdiv),
     $  h_r_wan(nrpts),
     $  tempmat1(num_orb_max,num_wancell),
     $  tempmat2(num_orb_max,num_orb_max),
     $  tempmat3(num_wancell,num_wancell),
     $  zinv(num_wancell,num_wancell),
     $  zinv_half(num_wancell,num_wancell),
     $  g_path_tmp(num_wancell,num_wancell,n_omega),
     $  g_path_tmp2,
     $  g_lat_kw(num_wancell,num_wancell),
     $  g_lat_inv_temp(num_wancell,num_wancell),      

     $  mat_loc_ackw(num_orb_max,num_orb_max),
     $  tempmat5(num_orb_max,num_orb_max),
     $  tempmat11(num_wancell, num_wancell),
     $  tempmat12(num_wancell, num_wancell),
     $  fac

      integer :: irequest,istatus(mpi_status_size)      
      logical :: mpiflag
      logical :: file_exist      
      
      double precision, allocatable :: tempdble(:,:),
     $  tdos(:), pdos(:,:),
     $  eigenvalues_fi_all(:,:)
      complex*16, allocatable :: tempmat4(:,:),
     $  g_path_trace(:,:),
     $  g_path_orb(:,:,:), 
     $  g_path(:,:,:,:),
     $  g_lat_kw_diagonal_orb(:,:,:)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     1. zone_folding information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ai=dcmplx(0.0d0, 1.0d0)
      pi=datan2(1.0d0, 1.0d0)*4.0d0

c$$$      if (num_cell .ne. 1) then
c$$$        
c$$$        gvec_conv_shift=0
c$$$        
c$$$        write(179,*) 'gvec_conv_shift cal'
c$$$        do ik=1, ndim_k(me1)
c$$$          irk=n_mpi_k(me1)+ik
c$$$          call cal_gvec_conv_shift(kvec(1,irk),
c$$$     $      gvec_conv_shift(1,1,irk))          
c$$$        enddo    
c$$$        
c$$$        call mpi_allreduce_int(gvec_conv_shift,
c$$$     $    3*numk*num_cell,mpi_sum,mpi_comm_world)
c$$$        if (maswrk) then
c$$$          write(iun,*) 'gvec_conv_shift'
c$$$          do irk=1, numk
c$$$            do icell=1, num_cell
c$$$              write(iun,'(5i5)') irk,icell,
c$$$     $          gvec_conv_shift(:,icell,irk)
c$$$            enddo
c$$$          enddo
c$$$        endif
c$$$
c$$$        write(179,*) 'gvec_conv_shift_fi test'
c$$$        do ik=1, ndim_k(me1)
c$$$          irk=n_mpi_k(me1)+ik
c$$$          do icell=1, num_cell
c$$$            do jcell=1, num_cell
c$$$              phaseval_dcmplx=0.0d0
c$$$              do igshift=1, num_cell
c$$$                phaseval_dcmplx(igshift)=phaseval_dcmplx(igshift)
c$$$     $            +cdexp(ai*
c$$$     $            sum(
c$$$     $            matmul(
c$$$     $            recip_lattice_conv,
c$$$     $            kvec(:,irk)+gvec_conv_shift(:,igshift,irk)
c$$$     $            )
c$$$     $            *
c$$$     $            (lpoint(:,icell)-lpoint(:,jcell))
c$$$     $            )
c$$$     $            )
c$$$              enddo
c$$$              write(179,'(3i5,20f12.6)') ik,icell,jcell,phaseval_dcmplx,
c$$$     $          cdabs(phaseval_dcmplx)
c$$$c$$$  if (cdabs(phaseval_dcmplx) .gt. 1.0d-6) then
c$$$c$$$  write(179,'(a,f12.6)')'wrong non zero off diagonal term'
c$$$c$$$  $            ,cdabs(phaseval_dcmplx)
c$$$c$$$  write(iun,'(a,f12.6)')'wrong non zero off diagonal term'
c$$$c$$$  $            ,cdabs(phaseval_dcmplx)
c$$$c$$$  endif
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$      endif
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     2. h_k in fine grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
! The next bit does three things:
! A. compute h_k_wan_co(q) = v_matrix(q)^\dagger *  eigenvalues(q) * v_matrix(q)
! B. Fourier interpolate from the q grid to the k grid.
! C. Distribute the Fourier-interpolated hamiltonian over nodes.
!   The result is in h_k_wan_fi, which is a num_wann x num_wann matrix.
! D. Diagonalize h_k_wan_fi.  The eigenvectors of h_k_wan_fi are in w_matrix_fi(distributed over nodes), and the eigenvalues are in eigenvalues_fi_all (the same on every node).
! h_k_wan_fi is the wannier90 hamiltonian interpolated onto the fine (k) grid.
! w_matrix_fi is the eigenvectors of the wannier90 hamiltonian interpolated onto the fine (k) grid.
! eigenvalues_fi_all is the eigenvalues of the wannier90 hamiltonian interpolated onto the fine (k) grid.

! num_wancell = num_wann

      allocate(h_k_wan_fi(num_wancell,num_wancell,ndim_k_max))
      h_k_wan_fi=0.0d0
! h_k_wan_fi is distributed over nodes.

      do iwan=1, num_wann
        do jwan=1, num_wann
          
          h_k_wan_co=0.0d0
          
          ! h_k_wan_co(q) = v_matrix(q)^\dagger *  eigenvalues(q) * v_matrix(q). This involves a sum over num_bands.
          do irk=1, nqdiv
            do ibnd=1, num_bands                
              h_k_wan_co(irk)=h_k_wan_co(irk)+
     $          eigenvalues(ibnd,irk)*v_matrix(ibnd,jwan,irk)
     $          *dconjg(v_matrix(ibnd,iwan,irk))
            enddo
          enddo
          
!     fourier transform to r space
        
          ! h_r_wan(r) = nqdiv^{-1} \sum_q h_r_wan(q) exp(-i 2 pi r.q)
          h_r_wan=0.0d0          
          do i1=1, nrpts
            do j1=1, nqdiv
              h_r_wan(i1)=h_r_wan(i1)
     $          +h_k_wan_co(j1)*
     $          cdexp(-2.0*ai*pi*
     $          sum(irvec(:,i1)*kpt_latt(:,j1)))
     $          /dble(nqdiv)
            enddo
          enddo

          if (num_cell .eq. 1) then !  num_cell=1

! h_k_wan_fi(k) = \sum_r h_r_wan(r) exp(i 2 pi k.r)
            do ik=1, ndim_k(me1)
              irk=n_mpi_k(me1)+ik
              do j1=1, nrpts
                h_k_wan_fi(iwan,jwan,ik)
     $            =h_k_wan_fi(iwan,jwan,ik)
     $            +h_r_wan(j1)*
     $            cdexp(2.0*ai*pi*
     $            sum(
     $            kvec(:,irk)
     $            *irvec(:,j1))
     $            )/dble(ndegen(j1))
              enddo
            enddo

c$$$          else
c$$$            do ik=1, ndim_k(me1)
c$$$              irk=n_mpi_k(me1)+ik
c$$$              do icell=1, num_cell
c$$$                do jcell=1, num_cell
c$$$                  iwancell=iwan+(icell-1)*num_wann
c$$$                  jwancell=jwan+(jcell-1)*num_wann
c$$$                  do igshift=1,num_cell ! gvec_conv_shift_fi vector
c$$$                    do j1=1, nrpts
c$$$                      phaseval
c$$$     $                  =sum(irvec(:,j1)
c$$$     $                  *matmul(
c$$$     $                  matmul(transpose(real_lattice),
c$$$     $                  recip_lattice_conv),
c$$$     $                  (kvec(:,irk)+gvec_conv_shift(:,igshift,irk))
c$$$     $                  )
c$$$     $                  )
c$$$                      phasetau
c$$$     $                  =sum((lpoint(:,icell)-lpoint(:,jcell))
c$$$     $                  *matmul(recip_lattice_conv,
c$$$     $                  (kvec(:,irk)+gvec_conv_shift(:,igshift,irk))
c$$$     $                  )
c$$$     $                  )
c$$$                      
c$$$                      h_k_wan_fi(iwancell,jwancell,ik)
c$$$     $                  =h_k_wan_fi(iwancell,jwancell,ik)
c$$$     $                  +h_r_wan(j1)*
c$$$     $                  cdexp(ai*(phaseval+phasetau))
c$$$     $                  /dble(ndegen(j1))/dble(num_cell)
c$$$                    enddo
c$$$                  enddo      
c$$$                enddo
c$$$              enddo    
c$$$            enddo
          endif !  (num_cell .eq. 1) then !  num_cell=1

        enddo ! jwan=1, num_wann
      enddo ! iwan=1, num_wann
  

      allocate(eigenvalues_fi(num_wancell,ndim_k_max)) ! num_wancell = num_wann
      eigenvalues_fi=0.0d0
      allocate(eigenvalues_fi_all(num_wancell,numk)) ! eigenvalues_fi_all is the same on all nodes.
      eigenvalues_fi_all=0.0d0
      allocate(w_matrix_fi(num_wancell,num_wancell,ndim_k_max))
      w_matrix_fi=0.0d0
      

! Calculate w_matrix_fi and eigenvalues_fi_all which are the eigenvectors and eigenvalues of h_k_wan_fi.
      w_matrix_fi=h_k_wan_fi

      eigenvalues_fi_all=0.0d0
      do ik=1, ndim_k(me1)
        irk=n_mpi_k(me1)+ik
        call hermitianeigen_cmplxdouble
     $    (num_wancell,eigenvalues_fi(1,ik),w_matrix_fi(1,1,ik))
        eigenvalues_fi_all(:,irk)=eigenvalues_fi(:,ik)
      enddo
      
      call mpi_allreduce_dble(eigenvalues_fi_all,
     $  numk*num_wancell,mpi_sum,mpi_comm_world)                    
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     renormalize projector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      allocate(proj
     $  (num_orb_max,num_wancell,num_atom,num_cell,ndim_k_max))
      allocate(proj_renorm
     $  (num_orb_max,num_wancell,num_atom,num_cell,ndim_k_max))

      
      proj=0.0d0
      proj_renorm=0.0d0      


! dmft_projector changes the indexing of the input w_matrix_fi, which is num_wann x num_wann x ndim_k_max.  The output matrix proj is num_orb_max x num_wann x num_atom x num_cell=1 x ndim_k_max.  The first index of w_matrix_fi is transformed into the first and third indices of proj.
! dmft_projector also applies a projector P_c on the left side of w_matrix_fi; the input matrix w_matrix_fi varies over all wanniers, and the output matrix proj contains only correlated wanniers.
! numk and kvec are not used.
! proj is the eigenvectors of the wannier90 hamiltonian interpolated onto the fine (k) grid, but with the left index turned into the first and third index.
! input: w_matrix_fi
! output: proj = P_c w_matrix_fi
      call dmft_projector
     $  (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,
     $  w_matrix_fi,proj)


      ! is_kpath means we are calculating A(k,omega)
      if (is_kpath .eq. 0) then

! cal_hmat_rhomat_from_proj calculates and prints out a hamiltonian
!   and a density matrix. It does not return anything to its caller.
        call cal_hmat_rhomat_from_proj
     $    (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,
     $    eigenvalues_fi,proj,
     $    ', bare projector, fine grid')

      endif

      allocate(nbnd_in_win(numk))
      allocate(bnd_in_win(num_wancell,numk))      
      
! find_bnd_in_win creates bnd_in_win which says which bands are inside the energy window [ewin_min,ewin_max] for at least one k-point, and nbnd_in_win which is the count of bands inside the window.  Even though bnd_in_win and nbnd_in_win are indexed over k-points, they have the same value for all k-points.
      call find_bnd_in_win(eigenvalues_fi_all,numk,num_wancell,
     $  ewin_min,ewin_max,den_mode,is_spinorbit,nbnd_in_win,bnd_in_win)
c$$$      if (maswrk) then
c$$$        write(iun,*) 'num_band information in a find grid'
c$$$        write(iun,'(a, 10000i5)') 'bands in the window',
c$$$     $    maxval(nbnd_in_win), nbnd_in_win
c$$$        do irk=1, numk
c$$$          do ibnd=1, nbnd_in_win(irk)
c$$$            jbnd=bnd_in_win(ibnd,irk)
c$$$            write(iun,'(2i5, 2f12.6)')
c$$$     $        ibnd, jbnd,
c$$$     $        eigenvalues_fi_all(jbnd,irk)
c$$$          enddo
c$$$        enddo
c$$$      endif
      
! projector_renormalization_allk returns proj_renorm =  ( ( \sqrt{proj * P_E *  proj^\dagger})^{-1/2} ) * proj * P_E.  proj_renorm has the property that proj_renorm proj_renorm^\dagger = 1.
! P_E is a projector which keeps only the sector defined by ind_bnd_k, which says which bands are inside an energy window.
! input: proj(num_orb_max,num_wann,num_atom,num_cell,ndim_kk_co_max). Structurally this means proj(orb+atom,wann,k)
!          proj = P_c w_matrix_fi
! output: proj_renorm(num_orb_max,num_wann,num_atom,num_cell,ndim_kk_co_max) Structurally this means proj_renorm(orb+atom,wann,k)
      call projector_renormalization_allk
     $  (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,    
     $  eigenvalues_fi,proj,nbnd_in_win,bnd_in_win,proj_renorm)


      ! is_kpath means we are calculating A((k,omega)
      if (is_kpath .eq. 0) then

! cal_hmat_rhomat_from_proj calculates and prints out a hamiltonian
!   and a density matrix. It does not return anything to its caller.
        call cal_hmat_rhomat_from_proj
     $    (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,
     $    eigenvalues_fi,proj_renorm,
     $    ', renormalized projector, fine grid')

! projector_compare does nothing unless (is_kpath .eq. 1), which is not true here.
! projector_compare calculates proj_orig * proj_ren^\dagger and prints the result.
        call projector_compare
     $    (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,          
     $    proj,proj_renorm,'fine_grid')

      endif

! proj is never used again after this point

! trans_basis.dat is created by comdmft.py's generate_initial_transformation method, unless a user supplies another and sets the appropriate option in comdmft.py. It contains a unitary matrix for each impurity - for instance a d-shell impurity would have a 5x5 unitary matrix.  comdmft.py sets trans_basis equal to the identity, so it has no effect on the calculation.
! In summary, trans_basis is equal to the identity and has no effect on the calculation.
      inquire(file="trans_basis.dat", exist=is_trans_basis)      
      
      if (is_trans_basis) then      
        allocate(trans_basis
     $    (num_orb_max,num_orb_max,n_iac_mat_i:n_iac_mat_f))
        trans_basis=0.0d0
        
! read trans_basis.dat, if it exists, into trans_basis
        if (maswrk) then
          open(unit=305, file='trans_basis.dat')
          allocate(tempdble(2,num_orb_max))
          do ii=n_iac_mat_i, n_iac_mat_f
            tempdble=0.0d0            
            if (ii .eq. 0) cycle
            do jorb=1, norb_iac(ii)
              read(305,*) ((tempdble(jj,iorb), jj=1, 2), 
     $          iorb=1,norb_iac(ii))
              do iorb=1, norb_iac(ii)
                trans_basis(jorb,iorb,ii)
     $            =dcmplx(tempdble(1,iorb), tempdble(2,iorb))
              enddo
            enddo
c$$$  test if it is unitary
            call check_if_unitary(norb_iac(ii),num_orb_max,
     $        trans_basis(1,1,ii),fail_ind)
            if (fail_ind .ne. 0) then
              write(iun,*) 'trans_basis matrix is not unitary'
              call ending                
            endif
          enddo
          deallocate(tempdble)
          close(305)
        endif
        call mpi_bcast(trans_basis,n_iac_mat*num_orb_max**2,
     $    mpi_double_complex,0,MPI_COMM_WORLD,mpierr)
      endif


      if (is_trans_basis) then
! In summary, trans_basis is equal to the identity and has no effect on the calculation.

! If trans_basis.dat exists, use it to rearrange the indexing of proj_renorm.
        do icell=1, num_cell !  ! num_cell=1
          do iatom=1, num_atom
            do ik=1, ndim_k(me1)
              irk=n_mpi_k(me1)+ik                          
              ii=iac(iatom,icell)
              tempmat1=proj_renorm(:,:,iatom,icell,ik)
              proj_renorm(:,:,iatom,icell,ik)=0.0d0              
              call zgemm('n','n',num_orb(iatom),num_wancell,
     $          num_orb(iatom),(1.d0,0.d0),
     $          trans_basis(1,1,ii),num_orb_max,
     $          tempmat1,num_orb_max,
     $          (0.d0,0.d0),
     $          proj_renorm(1,1,iatom,icell,ik),num_orb_max)
            enddo
          enddo
        enddo
      endif

! Done using trans_basis.dat.

! Done calculating proj_renorm; it will remain static after this.
! ! proj_renorm is a renormalized version of the eigenvectors of the wannier90 hamiltonian interpolated onto the fine (k) grid, but with the left index turned into the first and third index, and projected onto the correlated wannier orbitals via P_c. Possibly the indices have been further rearranged by multiplying by trans_basis.dat. The renormalization guarantees that proj_renorm proj_renorm^\dagger = 1.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     5. projected green's function at each frequency
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! First we allocate memory.
       ! is_kpath means we are calculating A((k,omega)
      if (is_kpath .eq. 0) then

        allocate(gloc(num_orb_max,num_orb_max,
     $    n_omega,n_iac_mat_i:n_iac_mat_f))
        allocate(projected_g(num_orb_max,num_orb_max,
     $    num_atom,num_cell,n_omega)) !
        allocate(delta_mat(num_orb_max,num_orb_max,n_omega,
     $    n_iac_mat_i:n_iac_mat_f))        
        allocate(e_projected_mat
     $    (num_orb_max,num_orb_max,n_iac_mat_i:n_iac_mat_f))
        gloc=0.0d0
        projected_g=0.0d0
        delta_mat=0.0d0
        e_projected_mat=0.0d0

        if (cal_mode .eq. 1) then ! doing a non-postprocessing run

          allocate(g_lat_kw_diagonal_orb
     $      (num_wancell,ndim_k_max,n_omega))
          g_lat_kw_diagonal_orb=0.0d0

        elseif (cal_mode .eq. 2) then ! calculating DOS after analytic continuation

          allocate(tdos(n_omega))
          allocate(pdos(n_omega, num_wancell))
          tdos=0.0d0
          pdos=0.0d0

        endif ! (cal_mode .eq. 2)

      else ! not (is_kpath .eq. 0)

        if (cal_mode .eq. 3) then ! calculating A(k,omega) after analytic continuation

          allocate(g_path_trace(numk,n_omega))
          allocate(g_path_orb(numk,n_omega,num_wancell))
          allocate(g_path(num_wancell,num_wancell,n_omega,ndim_k_max))
          g_path_trace=0.0d0
          g_path_orb=0.0d0          
          g_path=0.0d0
        endif

      endif ! not (is_kpath .eq. 0)
      
      allocate(h_nonint(num_wancell,num_wancell,ndim_k_max))
      h_nonint=0.0d0
! Done allocating memory.

! The point of the next block is to calculate h_nonint, and it would be trivial if zinv_m1_exist were false, i.e. if zinv_m1.dat did not exist, which is the case when doing dft+dmft.
! h_nonint = zinv_half * eigenvalues_fi * zinv_half
! zinv_half = \sqrt(1 + \sum_atoms proj_renorm^\dagger * zinv_m1 * proj_renorm)
! If zinv_m1.dat doesn't exist that means that zinv_m1 = 0 and h_nonint = eigenvalues_fi.
! This is the origin of zinv_m1.dat: run_dc runs ComDC, and then takes data from ComDC's sig_mat.dat and puts it zinv_m1_mat.  After ComLowH is run, delta_postprocessing/cal_zinv_m1_diagonal takes info from zinv_m1_mat.dat and puts it in zinv_m1.dat
      do ik=1, ndim_k(me1)

        irk=n_mpi_k(me1)+ik


        if (zinv_m1_exist) then

          ! zinv = 1
          zinv=0.0d0
          do iwancell=1, num_wancell
            zinv(iwancell,iwancell)=1.0d0
          enddo
          
           ! mat_loc_ackw = zinv_m1
          do icell=1, num_cell !  ! num_cell=1
            do iatom=1, num_atom
              ii=iac(iatom,icell)
              mat_loc_ackw=0.0d0
              do iorb=1, num_orb(iatom)
                do jorb=1, num_orb(iatom)
                  iio=ioac(iorb,jorb,ii)
                  if (iio .ne. 0) then
                    mat_loc_ackw(iorb,jorb)
     $                =mat_loc_ackw(iorb,jorb)
     $                +zinv_m1(iio,ii)
                  endif
                enddo
              enddo
              
! tempmat1 = zinv_m1 * proj_renorm
              tempmat1=0.0d0
              call zgemm('n','n',num_orb(iatom),num_wancell,
     $          num_orb(iatom),(1.d0,0.d0),
     $          mat_loc_ackw,num_orb_max,
     $          proj_renorm(1,1,iatom,icell,ik),num_orb_max,
     $          (0.d0,0.d0),
     $          tempmat1,num_orb_max)
              

! zinv = 1 + \sum_atoms proj_renorm^\dagger * zinv_m1 * proj_renorm
              call zgemm('c','n',num_wancell,num_wancell,
     $          num_orb(iatom),(1.d0,0.d0),
     $          proj_renorm(1,1,iatom,icell,ik),num_orb_max,
     $          tempmat1,num_orb_max,(1.d0,0.d0),
     $          zinv,num_wancell)
              
            enddo !  iatom=1, num_atom
          enddo ! icell=1, num_cell !  ! num_cell=1
          
          call hermitianeigen_cmplxdouble
     $      (num_wancell,diag,zinv)
! zinv now contains the eigenvectors of zinv, and diag contains the eigenvalues of zinv.

! zinv_half = \sqrt(zinv)
          zinv_half=0.0d0
          do iwancell=1, num_wancell ! num_wancell = num_wann
            do jwancell=1, num_wancell
              do kwancell=1, num_wancell                              
                zinv_half(iwancell,jwancell)
     $            =zinv_half(iwancell,jwancell)
     $            +zinv(iwancell,kwancell)*dsqrt(diag(kwancell))
     $            *dconjg(zinv(jwancell,kwancell))
              enddo
            enddo
          enddo
          
! h_nonint = zinv_half * eigenvalues_fi * zinv_half
          do iwancell=1, num_wancell
            do jwancell=1, num_wancell
              do kwancell=1, num_wancell                              
                h_nonint(iwancell,jwancell,ik)
     $            =h_nonint(iwancell,jwancell,ik)
     $            +zinv_half(iwancell,kwancell)
     $            *eigenvalues_fi(kwancell,ik)
     $            *zinv_half(kwancell,jwancell)
              enddo
            enddo
          enddo

        else ! not (zinv_m1_exist)

          do iwancell=1, num_wancell
            h_nonint(iwancell,iwancell,ik)=eigenvalues_fi(iwancell,ik)
          enddo

        endif ! not (zinv_m1_exist)

      enddo ! ik=1, ndim_k(me1)
      ! Done obtaining h_nonint, possibly with the help of zinv_m1.dat. Done with zinv_m1.

! This next block calculates g_lat_inv.
! g_lat_inv = frequency(iomega) - h_nonint - \sum_atoms proj_renorm^\dagger (sig(iomega) - dc - I * gamma_selfenergy) proj_renorm
! ----------------------
! This is the origin of dc.dat: If doing qsgw+dmft, then run_dc runs ComDC, and then takes data (dc[0,jj,kk]) from ComDC's sig_mat.dat and puts it dc_mat.dat.  Alternatively, if doing dft+dmft, then run_dc/cal_nominal_dc creates dc_mat.dat.  Alternatively, if ('dc_mat_to_read' in control), then run_dc only copies an old version of dc_mat.dat from somewhere and uses that. After ComLowH is run, delta_postprocessing/cal_dc_diagonal takes takes from dc_mat.dat and puts it in dc.dat.
! when doing lqsgw+dmft, during the first call to comlowh + comdc (before the loop starts), dc.dat is all zeros.  This is set up by comdmft.py's prepare_seed_dc_sig_and_wannier_dat.
! dc.dat is used both by ComLowH and by CTQMC.
! -----------------------
! Here is is the origin of sig = sig.dat:
! comdmft.py's run_impurity_solver runs CTQMC and EVALSIM, and then:
! reads in ctqmc's output from params.obs.json and updates convergence.log
! "self-energy" from params.obs.json is saved in sig_bare.dat
! "self-energy" is smoothed using Gaussian broadening and stored in sigma.
!    sigma is saved in  sig_smth.dat.
! If any element of imag(sigma) is positive, then sig_causality is set
!       to 0=False; otherwise it is 1=True.
! If any element of imag(sigma) is positive, then sigma_to_delta = sigma_old [read in from sig.dat].
!    If it is true then sigma_to_delta is a mix of sigma with sigma_old [ read in from sig.dat].
! sigma_to_delta is saved in sig.dat.
! In summary, sig.dat comes from params.obs.json, then gaussian smoothing, then a check for negative imaginary parts, then mixing.
! ------------
! On the initial self-consistent cycle of dft+dmft and lqsgw+dmft, things are different:
! comdmft.py's generate_initial_self_energy creates sig.dat.
! If ('initial_self_energy' in control,  comdmft.py copies that self_energy to sig.dat.
!   Possibly it also copies information from initial_impurity_dir.
! Otherwise, comdmft.py copies data from dc.dat to sig.dat.
! ----------------
! If doing lqsgw+dmft, then before doing the self-consistent cycle there is a comlowh+comdc run.  To prepare for this dc.dat is all zeros, and sig.dat contains all zero's and omega's.  This is set up by comdmft.py's prepare_seed_dc_sig_and_wannier_dat.
! -------------

      allocate(g_lat_inv(num_wancell,num_wancell,ndim_k_max,n_omega))
      inv_numk=1.0d0/dble(numk)
      g_lat_inv=0.0d0
      
      do iomega=1,n_omega
        
        do ik=1, ndim_k(me1)
          irk=n_mpi_k(me1)+ik

          do ibnd=1, num_wancell
            g_lat_inv(ibnd,ibnd,ik,iomega)=frequency(iomega)
          enddo

          g_lat_inv(:,:,ik,iomega)=
     $      g_lat_inv(:,:,ik,iomega)-h_nonint(:,:,ik)
        
!     embeding self-energy term:E(Sigma-DC)
          
          do icell=1, num_cell !  ! num_cell=1
            do iatom=1, num_atom
              ii=iac(iatom,icell)

              mat_loc_ackw=0.0d0

              do iorb=1, num_orb(iatom)
                do jorb=1, num_orb(iatom)
                  iio=ioac(iorb,jorb,ii)

                  if (iio .ne. 0) then
                    mat_loc_ackw(iorb,jorb)
     $                =mat_loc_ackw(iorb,jorb)
     $                +sig(iio,ii,iomega)
     $                -dc(iio,ii)
                  endif
                  
                  if (iorb .eq. jorb) then
                    mat_loc_ackw(iorb,jorb)
     $                =mat_loc_ackw(iorb,jorb)                    
     $                -dcmplx(0.0d0, gamma_selfenergy) ! add broadening
                  endif

                enddo ! jorb=1, num_orb(iatom)
              enddo ! iorb=1, num_orb(iatom)
              
              tempmat1=0.0d0
              call zgemm('n','n',num_orb(iatom),num_wancell,
     $          num_orb(iatom),(1.d0,0.d0),
     $          mat_loc_ackw,num_orb_max,
     $          proj_renorm(1,1,iatom,icell,ik),num_orb_max,
     $          (0.d0,0.d0),
     $          tempmat1,num_orb_max)
              

              call zgemm('c','n',num_wancell,num_wancell,
     $          num_orb(iatom),(-1.d0,0.d0),
     $          proj_renorm(1,1,iatom,icell,ik),num_orb_max,
     $          tempmat1,num_orb_max,(1.d0,0.d0),
     $          g_lat_inv(1,1,ik,iomega),num_wancell)
              
            enddo !  iatom=1, num_atom
          enddo !  icell=1, num_cell !  ! num_cell=1
          
        enddo !  do ik=1, ndim_k(me1)
      enddo ! do iomega=1,n_omega
      
! Either calculate the Fermi energy = ef or read it in from ef.dat, depending on is_recal_ef.
! ef.dat is set to zero by comdmft.py/prepare_initial_ef before the first call to ComLowH in self-consistent dft+dmft or qsgw+dmft.
! If doing self-consistent dft+dmft or qsgw+dmft, AND if not(is_metal==0), ComLowH never reads in ef.dat; it only recomputes it and saves it out during each iteration.
! If (is_metal==0), or if doing post-processing, ComLowH reads in ef.dat and uses its value.
! ef.dat is not used by other parts of ComDMFT.
c$$$  endif      
      if (is_recal_ef  .eq. 1) then

! todo. search_fermi_energy_brent calculates the new value of the Fermi level ef. Uses the 'doping' parameter in comdmft.ini, which turns into the delta_charge variable read in from comlowh.ini.
! -------------
! This routine finds a Fermi energy ef (used in gmat) that solves the equation 0=cal_tos_analytical_tail_fun.
! cal_tos_analytical_tail_fun=Tr(nimp(b))-nelectron
! nimp(b) = -2 \beta^{-1} N_k^{-1} \sum_k \sum_{\omega} \exp(-\beta \omega) diag[ gmat(k,\omega) - \omega^{-1} - g2mat(k,\omega_{max}) / \omega^2 - g3mat(k,\omega_{max}) / \omega^3] - N_k^{-1} \sum_k (-1/2 + g2mat(k,\omega_{max}) \beta/4)
! g2mat(k,\omega_{max}) = (g(k,\omega_{max}) + g^\dagger(k,\omega_{max}))/2 * \omega_{max}^2
! g3mat(k,\omega_{max})  = (g(k,\omega_{max}) - g^\dagger(k,\omega_{max})  - 2/\omega_{max})/2 * \omega_{max}^3
! gmat(k,\omega) = g(k,\omega) = (ef + frequency(\omega) - h_nonint - \sum_atoms proj_renorm^\dagger (sig(\omega) - dc - I * gamma_selfenergy) proj_renorm)^{-1}
! ---------
! g_gw is set to zero by the caller
! if no spin-orbit, nelectron=dble(nint(nelectron_dble*2))/2.0d0    +delta_charge/2.0d0
! if spin-orbit, nelectron=dble(nint(nelectron_dble))+delta_charge
! nelectron_dble=N_q^{-1} \sum_{q,b,o, eigenvalues(b,k) < 0} |v_matrix(b,o,k)|^2 
        call search_fermi_energy_brent
     $    (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,0)

c$$$        write(179,*) 'Ef=', ef

! This is a file write.
        if (maswrk) then
          open(unit=10, file='ef.dat')
          write(10,'(f19.12)') ef            
          close(10)
        endif

      else ! not (is_recal_ef  .eq. 1)

! This is a file read.
        if (maswrk) then
          open(unit=10, file='ef.dat', status='old')
          read(10,*) ef
          close(10)
        endif

        call mpi_bcast(ef,1,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)

      endif ! not (is_recal_ef  .eq. 1)


c$$$  write(179,*) 'e_projected', e_projected
       ! is_kpath means we are calculating A((k,omega)
      if (is_kpath .eq. 0) then

! This calculates and writes e_projected_mat.dat
! e_projected_mat is also used to calculate delta_mat.
! e_projected_mat(atom) = Nk^{-1} (ndeg_iac(atom))^{-1} \sum_k proj_renorm(atom) * h_nonint * proj_renorm(atom)^\dagger - ef
! e_projected_mat.dat is used by CTQMC after the following transformations:
! e_projected_mat.dat        ->projected_eig.dat->e_imp.dat->e_imp->e_imp_key->CTQMC

        do ik=1, ndim_k(me1)
          do icell=1, num_cell !  ! num_cell=1
            do iatom=1, num_atom
              ii=iac(iatom,icell)

! ndeg_iac = number of impurities (atoms) that are equivalent to each other
              fac=1.0d0/dble(numk)/dble(ndeg_iac(ii))              
              tempmat1=0.0d0

              call zgemm('n','n',num_orb(iatom),num_wancell,
     $          num_wancell,(1.d0,0.d0),
     $          proj_renorm(1,1,iatom,icell,ik),num_orb_max,
     $          h_nonint(1,1,ik),num_wancell,
     $          (0.d0,0.d0),
     $          tempmat1,num_orb_max)
              
              call zgemm('n','c',num_orb(iatom),num_orb(iatom),
     $          num_wancell,fac,
     $          tempmat1,num_orb_max,proj_renorm(1,1,iatom,icell,ik),
     $          num_orb_max,(1.0d0,0.0d0),
     $          e_projected_mat(1,1,ii),num_orb_max)

            enddo ! iatom=1, num_atom
          enddo ! icell=1, num_cell !  ! num_cell=1
        enddo !  ik=1, ndim_k(me1)
        
        call mpi_allreduce_dcmplx(e_projected_mat,
     $    n_iac_mat*num_orb_max**2,mpi_sum,mpi_comm_world)
        
        do ii=n_iac_mat_i, n_iac_mat_f
          if (ii .eq. 0) cycle          
          do iorb=1, norb_iac(ii)
            e_projected_mat(iorb,iorb,ii)
     $        =e_projected_mat(iorb,iorb,ii)-ef
          enddo
        enddo
            
        if (maswrk) then
          open(unit=305, file='e_projected_mat.dat')
          do ii=n_iac_mat_i, n_iac_mat_f
            if (ii .eq. 0) cycle
            do iorb=1, norb_iac(ii)
              write(305, '(1000(2f12.6,2x))')
     $          (e_projected_mat(iorb,jorb,ii),
     $          jorb=1, norb_iac(ii))
            enddo
          enddo
          close(305)
        endif ! maswrk

      endif ! (is_kpath .eq. 0)

      if (cal_mode .eq. 1) then  ! doing a non-postprocessing run
        g_lat_kw_diagonal_orb=0.0d0
      endif


      do iomega=1,n_omega
        
        do ik=1, ndim_k(me1)
          irk=n_mpi_k(me1)+ik

          g_lat_kw=0.0d0
          g_lat_inv_temp=g_lat_inv(:,:,ik,iomega)
          
          do ibnd=1, num_wancell
            g_lat_inv_temp(ibnd,ibnd)
     $        =g_lat_inv_temp(ibnd,ibnd)+ef
          enddo

! g_lat_kw = (g_lat_inv + ef)^{-1}
          call dcmplx_matinv(g_lat_inv_temp, g_lat_kw,
     $      num_wancell,num_wancell)

          ! is_kpath means we are calculating A((k,omega)
          if (is_kpath .eq. 0) then

! projected_g(atom,omega) = Nk^{-1} \sum_k proj_renorm * g_lat_kw * proj_renorm^\dagger
            do icell=1, num_cell !  ! num_cell=1
              do iatom=1, num_atom
                
                tempmat1=0.0d0
                tempmat2=0.0d0
                

                call zgemm('n','n',num_orb(iatom),num_wancell,
     $            num_wancell,(1.d0,0.d0),
     $            proj_renorm(1,1,iatom,icell,ik),num_orb_max,
     $            g_lat_kw,num_wancell,(0.d0,0.d0),
     $            tempmat1,num_orb_max)
                
                call zgemm('n','c',num_orb(iatom),num_orb(iatom),
     $            num_wancell,(1.0d0,0.d0),
     $            tempmat1, num_orb_max,
     $            proj_renorm(1,1,iatom,icell,ik),
     $            num_orb_max,(0.0d0,0.d0),tempmat2,num_orb_max)
                
                do iorb=1, num_orb(iatom)
                  do jorb=1, num_orb(iatom)
                    projected_g(iorb,jorb,iatom,icell,iomega)
     $                =projected_g(iorb,jorb,iatom,icell,iomega)
     $                +tempmat2(iorb,jorb)
     $                /dble(numk)
                  enddo
                enddo

              enddo !  iatom=1, num_atom
            enddo !  icell=1, num_cell !  ! num_cell=1

! g_lat_kw_diagonal_orb(k,omega) = diag (w_matrix_fi * g_lat_kw * w_matrix_fi^\dagger)
            if (cal_mode .eq. 1) then ! doing a non-postprocessing run
              do iwancell=1, num_wancell
                do jwancell=1, num_wancell
                  do kwancell=1, num_wancell
                    g_lat_kw_diagonal_orb(iwancell,ik,iomega)
     $                =g_lat_kw_diagonal_orb(iwancell,ik,iomega)
     $                +w_matrix_fi(iwancell,jwancell,ik)
     $                *g_lat_kw(jwancell,kwancell)
     $                *dconjg(w_matrix_fi(iwancell,kwancell,ik))
                  enddo
                enddo
              enddo

            elseif (cal_mode .eq. 2) then ! calculating DOS after analytic continuation

! tdos(iomega) = - Nk^{-1} pi^{-1} \sum_k Tr(imag(g_lat_kw))
              do iwancell=1, num_wancell ! num_wancell = num_wann
                tdos(iomega)=tdos(iomega)
     $            -dimag(g_lat_kw(iwancell,iwancell))
     $            /pi/dble(numk)
                
! pdos(iomega) = -Nk^{-1} pi^{-1} \sum_k imag(w_matrix_fi * g_lat_kw * w_matrix_fi^\dagger)
                do jwancell=1, num_wancell
                  do kwancell=1, num_wancell
                    pdos(iomega,iwancell)=pdos(iomega,iwancell)
     $                -dimag(w_matrix_fi(iwancell,jwancell,ik)
     $                *g_lat_kw(jwancell,kwancell)
     $                *dconjg(w_matrix_fi(iwancell,kwancell,ik)))
     $                /pi/dble(numk)
                  enddo
                enddo

              enddo !  num_wancell = num_wann

            endif ! (cal_mode .eq. 2)

! g_path_orb(k,omega) = I * (2 pi)^{-1} diag(w_matrix_fi * (g_lat_kw - g_lat_kw^\dagger) w_matrix_fi^\dagger)
! g_path_trace(k,omega) = I * (2 pi)^{-1} trace(w_matrix_fi * (g_lat_kw - g_lat_kw^\dagger) w_matrix_fi^\dagger)
! The following elseif line refers to the  if(is_kpath .eq. 0)  statement; it executes if not(is_kpath .eq. 0).
          elseif (cal_mode .eq. 3) then ! calculating A(k,omega) after analytic continuation
            
            tempmat11=(g_lat_kw-transpose(dconjg(g_lat_kw)))*ai/2.0d0/pi
            call zgemm('n','n',num_wancell,num_wancell,
     $        num_wancell,(1.d0,0.d0),
     $        w_matrix_fi(1,1,ik),num_wancell,
     $        tempmat11,num_wancell,
     $        (0.d0,0.d0),
     $        tempmat12,num_wancell)
            
            call zgemm('n','c',num_wancell,num_wancell,
     $        num_wancell,(1.d0,0.d0),
     $        tempmat12,num_wancell,w_matrix_fi(1,1,ik),num_wancell,
     $        (1.d0,0.d0),
     $        g_path(1,1,iomega,ik),num_wancell)            
            
            do iwancell=1, num_wancell
              g_path_trace(irk,iomega)
     $          =g_path_trace(irk,iomega)
     $          +g_path(iwancell,iwancell,iomega,ik)
              
              g_path_orb(irk,iomega,iwancell)
     $          =g_path_orb(irk,iomega,iwancell)
     $          +g_path(iwancell,iwancell,iomega,ik)
            enddo

! The following endif line refers to the  if(is_kpath .eq. 0)  statement; it executes if not(is_kpath .eq. 0).
          endif

        enddo !  ik=1, ndim_k(me1)
      enddo ! iomega=1,n_omega
      

      ! is_kpath means we are calculating A((k,omega)
      if (is_kpath .eq. 0) then

        
        call mpi_allreduce_dcmplx(projected_g,
     $    n_omega*num_orb_max**2*num_atom*num_cell, !  ! num_cell=1
     $    mpi_sum,mpi_comm_world)

! calculate and write out g_loc_mat.dat
! gloc = projected_g / ndeg_iac
! # g_loc_mat.dat is converted to g_loc.dat and used by ComDC.
# It is not used if (control['dc_mode'] == 'dc_scf' not 'dc_at_gw') and (control['dc_g'] == 'gimp' not 'gloc').
# In this case CTQMC's gimp.dat is used instead.

! This copies projected_g into gloc, and divides by ndeg_iac.
        gloc=0.0d0
        do iomega=1, n_omega
          do icell=1,num_cell !  ! num_cell=1
            do iatom=1, num_atom
              do iorb=1, num_orb(iatom)
                do jorb=1, num_orb(iatom)
                  gloc(iorb,jorb,iomega,iac(iatom,icell))
     $              =gloc(iorb,jorb,iomega,iac(iatom,icell))
     $              +projected_g(iorb,jorb,iatom,icell,iomega)
     $              /dble(ndeg_iac(iac(iatom,icell)))
                enddo
              enddo
            enddo
          enddo
        enddo
        
        
        if (maswrk) then

c$$$          open(unit=305, file='g_loc.dat')
c$$$          do iomega=1, n_omega
c$$$            tempdcmplx=0.0d0
c$$$            write(305,'(f19.12, 3x)', advance='no') 
c$$$     $        frequency_dble(iomega)
c$$$            
c$$$            do icell=1,num_cell
c$$$              do iatom=1, num_atom
c$$$                ii=iac(iatom,icell)                            
c$$$                do iorb=1, norb_iac(ii)
c$$$                  do jorb=1, norb_iac(ii)
c$$$                    iio=ioac(iorb,jorb,ii)
c$$$                    if (iio .eq. 0) cycle                                  
c$$$                    tempdcmplx(iio,ii)
c$$$     $                =tempdcmplx(iio,ii)
c$$$     $                +projected_g(iorb,jorb,iatom,icell,iomega)
c$$$     $                /dble(ndeg_itot(iio,ii))                  
c$$$                  enddo
c$$$                enddo
c$$$              enddo
c$$$            enddo
c$$$            
c$$$            do ii=n_iac_mat_i, n_iac_mat_f
c$$$              if ((is_magnetic .eq. 1) .and. (ii .eq. 0)) cycle
c$$$              do iio=1, maxval(ioac(:,:,ii))
c$$$                write(305,'(2f19.12, 2x)', advance='no')
c$$$     $            tempdcmplx(iio,ii)
c$$$              enddo
c$$$            enddo
c$$$            write(305,*)
c$$$          enddo
c$$$          close(305)
          
          open(unit=305, file='g_loc_mat.dat')
          do iomega=1, n_omega
            write(305,'(f19.12, 3x)', advance='no') 
     $        frequency_dble(iomega)
            
            do ii=n_iac_mat_i, n_iac_mat_f
              if ((is_magnetic .eq. 1) .and. (ii .eq. 0)) cycle
              do jorb=1, norb_iac(ii)
                do iorb=1, norb_iac(ii)
                  write(305,'(2f19.12, 2x)', advance='no')
     $              gloc(iorb,jorb,iomega,ii)
                enddo
              enddo
            enddo
            write(305,*)
          enddo
          close(305)
          
        endif ! maswrk

! calculate and write out delta_mat.dat
! delta_mat(omega) = \sum_atoms (frequency(iomega) - e_projected_mat - (projected_g)^{-1} + dc - sig(omega))/ndeg(atom)
! delta_mat.dat is used by CTQMC.
! -------------
! This is the origin of dc.dat: If doing qsgw+dmft, then run_dc runs ComDC, and then takes data (dc[0,jj,kk]) from ComDC's sig_mat.dat and puts it dc_mat.dat.  Alternatively, if doing dft+dmft, then run_dc/cal_nominal_dc creates dc_mat.dat.  Alternatively, if ('dc_mat_to_read' in control), then run_dc only copies an old version of dc_mat.dat from somewhere and uses that. After ComLowH is run, delta_postprocessing/cal_dc_diagonal takes takes from dc_mat.dat and puts it in dc.dat.
! when doing lqsgw+dmft, during the first call to comlowh + comdc (before the loop starts), dc.dat is all zeros.  This is set up by comdmft.py's prepare_seed_dc_sig_and_wannier_dat.
! dc.dat is used both by ComLowH and by CTQMC.
! ------------
! Here is is the origin of sig=sig.dat:
! comdmft.py's run_impurity_solver runs CTQMC and EVALSIM, and then:
! reads in ctqmc's output from params.obs.json and updates convergence.log
! "self-energy" from params.obs.json is saved in sig_bare.dat
! "self-energy" is smoothed using Gaussian broadening and stored in sigma.
!    sigma is saved in  sig_smth.dat.
! If any element of imag(sigma) is positive, then sig_causality is set
!       to 0=False; otherwise it is 1=True.
! If any element of imag(sigma) is positive, then sigma_to_delta = sigma_old [read in from sig.dat].
!    If it is true then sigma_to_delta is a mix of sigma with sigma_old [ read in from sig.dat].
! sigma_to_delta is saved in sig.dat.
! In summary, sig.dat comes from params.obs.json, then gaussian smoothing, then a check for negative imaginary parts, then mixing.
! ------------
! On the initial self-consistent cycle of dft+dmft and lqsgw+dmft, things are different:
! comdmft.py's generate_initial_self_energy creates sig.dat.
! If ('initial_self_energy' in control,  comdmft.py copies that self_energy to sig.dat.
!   Possibly it also copies information from initial_impurity_dir.
! Otherwise, comdmft.py copies data from dc.dat to sig.dat.
! ----------------
! If doing lqsgw+dmft, then before doing the self-consistent cycle there is a comlowh+comdc run.  To prepare for this dc.dat is all zeros, and sig.dat contains all zero's and omega's.  This is set up by comdmft.py's prepare_seed_dc_sig_and_wannier_dat.
! -------------
        do iomega=1, n_omega
          
          do icell=1, num_cell ! num_cell=1
            do iatom=1, num_atom

! tempmat2 = (projected_g)^{-1}
              tempmat2=0.0d0
              call dcmplx_matinv
     $          (projected_g(1,1,iatom,icell,iomega),
     $          tempmat2, num_orb(iatom), num_orb_max)

              ii=iac(iatom,icell)
              if ((is_magnetic .eq. 1) .and. (ii .eq. 0)) cycle                            
              
              do iorb=1, norb_iac(ii)
                do jorb=1, norb_iac(ii)
                  iio=ioac(iorb,jorb,ii)

                  if (iorb .eq. jorb) then
                    deltaij=1.0d0
                  else
                    deltaij=0.0d0
                  endif

! delta_mat(omega) += \sum_atoms (frequency(iomega) - e_projected_mat - (projected_g)^{-1} + dc - sig(omega))/ndeg(atom)

                  delta_mat(iorb,jorb,iomega,ii)
     $              =delta_mat(iorb,jorb,iomega,ii)
     $              +(frequency(iomega)*deltaij
     $              -e_projected_mat(iorb,jorb,ii)
     $              -tempmat2(iorb,jorb))/ndeg_iac(ii)

                  if (iio .eq. 0) cycle

                  delta_mat(iorb,jorb,iomega,ii)
     $              =delta_mat(iorb,jorb,iomega,ii)
     $              +(dc(iio,ii)-sig(iio,ii,iomega))/ndeg_iac(ii)

                enddo ! jorb=1, norb_iac(ii)
              enddo ! iorb=1, norb_iac(ii)

            enddo ! iatom=1, num_atom
          enddo ! icell=1, num_cell ! num_cell=1
        enddo ! iomega=1, n_omega
        

        if (maswrk) then
          open(unit=305, file='delta_mat.dat')
          do iomega=1, n_omega
            write(305,'(f19.12, 3x)', advance='no') 
     $        frequency_dble(iomega)
            
            do ii=n_iac_mat_i, n_iac_mat_f
              if ((is_magnetic .eq. 1) .and. (ii .eq. 0)) cycle
              do jorb=1, norb_iac(ii)
                do iorb=1, norb_iac(ii)
                  write(305,'(2f19.12, 2x)', advance='no')
     $              delta_mat(iorb,jorb,iomega,ii)
                enddo
              enddo
            enddo
            write(305,*)
          enddo ! iomega=1, n_omega
          close(305)
        endif ! maswrk


        if (cal_mode .eq. 1) then ! doing a non-postprocessing run

! fermi_level_spectral_function_local calculates and prints out local_spectral_matrix_ef.dat, which is under certain circumstances used by comdmft.py. It also writes out n_loc.dat, which is never used by ComDMFT.
          call fermi_level_spectral_function_local

! ! fermi_level_spectral_function_ksum computes and prints out  spectral_matrix_ef.dat
          call fermi_level_spectral_function_ksum
     $      (numk, ndim_k, n_mpi_k, ndim_k_max, g_lat_kw_diagonal_orb)
          
          allocate(eigenvalues_co(num_wancell,ndim_kk_co_max)) !
          eigenvalues_co=0.0d0
          allocate(eigenvalues_co_all(num_wancell,nqdiv))
          eigenvalues_co_all=0.0d0
          allocate(eigenvalues_co_rest(num_rest,nqdiv))
          eigenvalues_co_rest=0.0d0          
          allocate(w_matrix_co(num_wancell,num_wancell,ndim_kk_co_max))
          w_matrix_co=0.0d0
          
          allocate(vw_matrix_co(num_bands,num_bands,nqdiv))
          vw_matrix_co=0.0d0          
          
          allocate(proj_co
     $      (num_orb_max,num_wancell,num_atom,num_cell,ndim_kk_co_max))
          allocate(proj_renorm_co
     $      (num_orb_max,num_wancell,num_atom,num_cell,ndim_kk_co_max))
          
          
          proj_co=0.0d0
          proj_renorm_co=0.0d0

          allocate(h_nonint_co(num_wancell,num_wancell,ndim_kk_co_max))
          h_nonint_co=0.0d0
          allocate(to_orig_bnd(num_bands,nqdiv))
          to_orig_bnd=0
          allocate(n_matrix(num_bands,num_bands,nqdiv))
          n_matrix=0.0d0
          
! cal_mat_eig_co calculates h_nonint_co, and eigenvalues_co, and eigenvalues_co_all, and w_matrix_co. It would be trivial if zinv_m1_exist were false, i.e. if zinv_m1.dat did not exist, which is the case when doing dft+dmft.
! h_nonint_co is used only in cal_n_matrix, and w_matrix_co is used by write_den_matrix.
! eigenvalues_co=eigenvalues_co_all and w_matrix_co are the eigenvalues and eigenvectors of v_matrix^\dagger * eigenvalues * v_matrix
! h_nonint_co(k) = \sqrt{zinv}(k) * eigenvalues_co(k) * \sqrt{zinv}(k)
! zinv = 1 + \sum_i proj_renorm_co^\dagger(i,k) * zinv_m1(i) * proj_renorm_co(i,k)
! If zinv_m1.dat doesn't exist that means that zinv_m1 = 0 and h_nonint_co = eigenvalues_co.
! output: h_nonint_co
! output: eigenvalues_co, eigenvalues_co_all
! output: w_matrix_co
! input: v_matrix
! input: eigenvalues
! input: zinv_m1
! input: proj_renorm_co
          call cal_mat_eig_co
   
! dmft_projector changes the indexing of the input w_matrix_co(num_wann,num_wann,k) to proj_co(orb+atom,num_wann,k).
! dmft_projector also applies a projector P_c on the left side of w_matrix_co; the input matrix w_matrix_co varies over all wanniers, and the output matrix proj_co contains only correlated wanniers.
! numk and kvec are not used.
! input: w_matrix_co
! output: proj_co = P_c w_matrix_co
          call dmft_projector
     $      (nqdiv,ndim_kk_co,n_mpi_kk_co,ndim_kk_co_max,kpt_latt,
     $      w_matrix_co,proj_co)

          allocate(nbnd_in_win_co(nqdiv))
          allocate(bnd_in_win_co(num_wancell,nqdiv))

! bnd_in_win_co says which bands are inside the energy window [ewin_min,ewin_max] for at least one q-point, and nbnd_in_win_co which the count of bands inside the window.
! In summary, bnd_in_win_co defines a projection operator P_E.
! The difference about the _co variables is that they are defined on nqdiv q-points, not numk k-points.
! den_mode = 'bnd'
          if (den_mode .eq. 'erg') then
            call find_bnd_in_win(eigenvalues_co_all,nqdiv,num_wancell,
     $        ewin_min,ewin_max,den_mode,is_spinorbit,
     $        nbnd_in_win_co,bnd_in_win_co)
          elseif (den_mode .eq. 'bnd') then
            do irk=1, nqdiv
              nbnd_in_win_co(irk)=nbnd_in_win(1)
              do ibnd=1, nbnd_in_win(1)
                bnd_in_win_co(ibnd,irk)=bnd_in_win(ibnd,1)
              enddo
            enddo
          endif
          
! projector_renormalization_allk returns proj_renorm_co =  ( ( \sqrt{proj_co * P_E *  proj_co^\dagger})^{-1/2} ) * proj_co * P_E.  proj_renorm_co has the property that proj_renorm_co proj_renorm_co^\dagger = 1.
! P_E is a projector which keeps only the sector defined by ind_bnd_k, which says which bands are inside an energy window.
! input: proj_co(num_orb_max,num_wann,num_atom,num_cell,ndim_kk_co_max). Structurally this means proj_co(orb+atom,wann,k)
!        proj_co = P_c w_matrix_co
! output: proj_renorm_co(num_orb_max,num_wann,num_atom,num_cell,ndim_kk_co_max) Structurally this means proj_renorm_co(orb+atom,wann,k)
          call projector_renormalization_allk
     $      (nqdiv,ndim_kk_co,n_mpi_kk_co,ndim_kk_co_max,kpt_latt,
     $      eigenvalues_co,proj_co,nbnd_in_win_co,bnd_in_win_co,
     $      proj_renorm_co)
          
! projector_compare does nothing unless (is_kpath .eq. 1), which is not true because cal_mode=1.
! projector_compare calculates proj_orig * proj_ren^\dagger and prints the result.
          call projector_compare
     $      (nqdiv,ndim_kk_co,n_mpi_kk_co,ndim_kk_co_max,kpt_latt,
     $      proj_co,proj_renorm_co,'coarse_grid')

! no further use of proj_co after this

          if (is_trans_basis) then
! In summary, trans_basis is equal to the identity and has no effect on the calculation.

            do icell=1, num_cell ! num_cell=1
              do iatom=1, num_atom
                do ik=1, ndim_kk_co(me1)
                  irk=n_mpi_kk_co(me1)+ik                          
                  ii=iac(iatom,icell)
                  tempmat1=proj_renorm_co(:,:,iatom,icell,ik)
                  proj_renorm_co(:,:,iatom,icell,ik)=0.0d0              
                  call zgemm('n','n',num_orb(iatom),num_wancell,
     $              num_orb(iatom),(1.d0,0.d0),
     $              trans_basis(1,1,ii),num_orb_max,
     $              tempmat1,num_orb_max,
     $              (0.d0,0.d0),
     $              proj_renorm_co(1,1,iatom,icell,ik),num_orb_max)
                enddo
              enddo
            enddo
          endif

! cal_vw_matrix_co calculates vw_matrix_co and to_orig_bnd, which are used only by write_den_matrix. It also calculates eigenvalues_co_rest, which is used only by cal_n_matrix.
          call cal_vw_matrix_co

! cal_n_matrix calculates the n_matrix:
! It calculates g_lat_inv_co(omega) = frequency(iomega) - h_nonint_co(:,:,ik) + \sum_{atom} proj_renorm_co(k)^\dagger * (sig(omega) - dc - I * gamma_selfenergy) * proj_renorm_co(k)
! proj_renorm_co is a renormalized version of proj_co = P_c w_matrix_co
! Next it calculates gmat = (ef + g_lat_inv_co(:,:,n_omega))^{-1}
! Next it uses gmat   to calculate g2mat and g3mat, which will help with the tail at large momenta.
! Next, it loops over omega, and in this loop it recalculates gmat = (ef + g_lat_inv_co(:,:,omega))^{-1}, and then uses gmat, g2mat, g3mat to add a contribution to n_matrix.
! All of this logic has updated only the P_wann sector of n_matrix. Next we update the P_rest sector of of n_matrix, by adding   fermi(beta(eigenvalues_co_rest(k) -ef)), where eigenvalues_co_rest is a diagonal matrix and is defined only in the P_rest sector.
! Last it repeats the above in a loop over k, so n_matrix is summed over k.
! ---------------
! output: n_matrix(num_bands,num_bands,nqdiv)
! local: g_lat_inv_co(num_wann,num_wann,n_omega)
! input: h_nonint_co(num_wann,num_wann,ndim_kk_co_max)
! input: proj_renorm_co(num_orb_max,num_wann,num_atom,num_cell,ndim_kk_co_max)
! proj_renorm_co is a renormalized version of proj_co = P_c w_matrix_co
! input: eigenvalues_co_rest(num_rest,nqdiv)
          call cal_n_matrix

! no further use of proj_renorm_co

! write_den_matrix writes out wannier_den_matrix.dat, which is used by rspflapw.
! wannier_den_matrix.dat can be in either binary or hdf5 format.
! It writes:
!       - "temperature"
!       - "nqdiv" - the number of q-points
!       - "n_dmft_bnd" from nbnd_in_win_co [called 'num_band_in' when saved from comrisb, and read in by rspflapw as 'n_dmft_bnd'] - should be the number of bands
!       - "n_dmft_wan" from nbnd_in_win_co
!       - "include_band_new"
!               include_band_new comes from nbnd_in_win_co, bnd_in_win_co, to_orig_bnd, and include_bands.
!       - "vw_matrix_new"
!               vwmatrix_new is a subsector of vw_matrix_co, with the sector specified by nbnd_in_win_co, to_orig_bnd, bnd_in_win_co. vwmatrix_new is renormalized so that it satisfies vwmatrix_new^\dagger vwmatrix_new = 1.
!       - 'n_matrix_new' - the density matrix
!                n_matrix_new equal is a subsector of n_matrix, with the sector specified by nbnd_in_win_co, bnd_in_win_co.
          call write_den_matrix


          deallocate(n_matrix)
          deallocate(vw_matrix_co)
          deallocate(to_orig_bnd)
          deallocate(nbnd_in_win_co)
          deallocate(bnd_in_win_co)          
          deallocate(h_nonint_co)
          
          deallocate(eigenvalues_co)
          deallocate(eigenvalues_co_all)
          deallocate(eigenvalues_co_rest)                    
          deallocate(w_matrix_co)
          deallocate(proj_co)
          deallocate(proj_renorm_co)


        elseif (cal_mode .eq. 2) then ! calculating DOS after analytic continuation
          call mpi_allreduce_dble(tdos,
     $      n_omega,mpi_sum,mpi_comm_world)
          call mpi_allreduce_dble(pdos,
     $      n_omega*num_wancell,mpi_sum,mpi_comm_world)            
          if (maswrk) then          
            open(unit=305, file='tdos.dat')
            do iomega=1,n_omega
              write(305,'(2f19.12)') 
     $          frequency_dble(iomega),tdos(iomega)
            enddo
            close(305)
            open(unit=305, file='pdos.dat')
            do iomega=1,n_omega
              write(305,'(f19.12, 5x, 300(f19.12,2x))') 
     $          frequency_dble(iomega),(pdos(iomega, iwancell),
     $          iwancell=1, num_wancell)
            enddo
            close(305)
          endif
          
c$$$  if (num_cell .eq. 1) then
          inquire(file="momentum_optics_wan.dat", exist=file_exist)
          if (file_exist) then
            if (maswrk) call read_wann_inip
! This is the only place where the inip file is read by comlowh, and the only place where it is used.
            call mpi_bcast(num_alimj,1,
     $        mpi_integer,0,MPI_COMM_WORLD,mpierr)!
            if (.not. maswrk) allocate(nalimj(num_wann))
            call mpi_bcast(nalimj,num_wann,
     $        mpi_integer,0,MPI_COMM_WORLD,mpierr)!
            call mpi_bcast(max_nalimj,1,
     $        mpi_integer,0,MPI_COMM_WORLD,mpierr)!
            if (.not.maswrk)
     $        allocate(alimj(num_alimj,max_nalimj,num_wann))
            call mpi_bcast(alimj,num_alimj*max_nalimj*num_wann,
     $        mpi_integer,0,mpi_comm_world,mpierr) !           
            if (.not.maswrk)
     $        allocate(alimj_coeff(max_nalimj,num_wann))
            call mpi_bcast(alimj_coeff,max_nalimj*num_wann,
     $        mpi_double_complex,0,mpi_comm_world,mpierr)!            
            if (.not.maswrk)
     $        allocate(alimj_xaxis(3,max_nalimj,num_wann))
            call mpi_bcast(alimj_xaxis,3*max_nalimj*num_wann,
     $        mpi_double_precision,0,mpi_comm_world,mpierr)!
            if (.not.maswrk)
     $        allocate(alimj_zaxis(3,max_nalimj,num_wann))
            call mpi_bcast(alimj_zaxis,3*max_nalimj*num_wann,
     $        mpi_double_precision,0,mpi_comm_world,mpierr) !           
c$$$            write(179,*) num_alimj
c$$$            write(179,*) nalimj
c$$$            write(179,*) max_nalimj
c$$$            write(179,*) alimj_coeff
c$$$            write(179,*) alimj_xaxis
c$$$            write(179,*) alimj_zaxis
            
            call prepare_optics

            deallocate(nalimj)
            deallocate(alimj)
            deallocate(alimj_coeff)
            deallocate(alimj_xaxis)
            deallocate(alimj_zaxis)
            endif ! (file_exist)
c$$$  endif
          endif !  (cal_mode .eq. 2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     5. impurity level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This elseif(cal_mode .eq. 3) continues a long-ago if (is_kpath .eq. 0) which covers cal_mode=1,2
      elseif (cal_mode .eq. 3) then ! calculating A(k,omega) after analytic continuation
        call mpi_allreduce_dcmplx(g_path_trace,
     $    n_omega*numk,mpi_sum,mpi_comm_world)
        call mpi_allreduce_dcmplx(g_path_orb,
     $    n_omega*numk*num_wancell,mpi_sum,mpi_comm_world)        
        if (maswrk) then
          open(unit=990, file='spectral.matrix',form = 'unformatted')
          write(990) is_spinorbit,1,numk,num_wancell
          write(990) n_omega
          write(990) frequency_dble
        endif
        
        if (.not. maswrk) then
          do ik=1, ndim_k(me1)
            irk=n_mpi_k(me1)+ik
            call mpi_isend(g_path(1,1,1,ik),
     $        num_wancell**2*n_omega,mpi_double_complex,0,irk,
     $        mpi_comm_world,irequest,mpierr)
            call mpi_wait(irequest,istatus,mpierr)
            call mpi_test(irequest,mpiflag,istatus,mpierr)
            if (.not.mpiflag)
     $        call mpi_abort(mpi_comm_world,mpierr)
          enddo
        else ! maswrk
          do irk=1, numk
             find_proc_ind(ndim_k, n_mpi_k, nproc,irk,
     $        proc_num,indout)
            if (proc_num .eq. 0)then
              write(990) g_path(:,:,:,indout)
            else
              g_path_tmp=0.0d0
              call mpi_irecv(g_path_tmp,
     $          num_wancell**2*n_omega,mpi_double_complex,proc_num,
     $          irk,mpi_comm_world,irequest,mpierr)
              call mpi_wait(irequest,istatus,mpierr)
              call mpi_test(irequest,mpiflag,istatus,mpierr)
              if (.not.mpiflag)
     $          call mpi_abort(mpi_comm_world,mpierr)
              write(990) g_path_tmp
              do iomega=1,n_omega
                g_path_tmp2=0.0d0
                do ibnd=1, num_wancell
                  g_path_tmp2
     $              =g_path_tmp2
     $              +g_path_tmp(ibnd,ibnd,iomega)
                enddo
                if (
     $            (dble(g_path_tmp2)
     $            -dble(g_path_trace(irk,iomega)))
     $            .gt. 1.0d-6)  then
                  write(iun,'(a,2i6,2f20.12)')
     $              'something wrong in communication',
     $              irk,iomega,dble(g_path_tmp2),
     $              dble(g_path_trace(irk,iomega))
                endif
              enddo
            endif
          enddo
        endif ! maswrk

        call mpi_barrier(mpi_comm_world,mpierr)
        if (maswrk) then
          close(990)
        endif
        
        if (maswrk) then          
          open(unit = 990, file = 'spectral.dat',
     $      form = 'formatted')
          do iomega=1, n_omega
            do irk=1, numk
              write(990,'(i7, 2(f20.12,2x))')
     $          irk, frequency_dble(iomega),
     $          dble(g_path_trace(irk,iomega))
            enddo
          enddo
          close(990)

          open(unit = 990, file = 'spectral_orb.dat',
     $      form = 'formatted')
          do iomega=1, n_omega
            do irk=1, numk
              write(990,'(i7, (f20.12,2x),1000(f20.12,2x))')
     $          irk, frequency_dble(iomega),
     $          (dble(g_path_orb(irk,iomega,iwancell)),
     $          iwancell=1, num_wancell)
            enddo
          enddo
          close(990)

          open(unit = 990, file = 'wannier_band_non_interpolated.dat',
     $      form = 'formatted')
          do ibnd=1,num_wancell
            do irk=1, numk
              write(990,'(2i7, f20.12)')
     $          irk, ibnd,eigenvalues_fi_all(ibnd,irk)
            enddo
          enddo
          close(990)
        endif ! maswrk
      endif ! (cal_mode .eq. 3)


      deallocate(nbnd_in_win)
      deallocate(bnd_in_win)
      
      deallocate(h_k_wan_fi)
      deallocate(eigenvalues_fi)
      deallocate(eigenvalues_fi_all)      
      deallocate(w_matrix_fi)
      deallocate(h_nonint)
      deallocate(proj)
      deallocate(proj_renorm)
      deallocate(g_lat_inv)
      if (is_trans_basis) then            
        deallocate(trans_basis)
      endif
      
      ! is_kpath means we are calculating A((k,omega)
      if (is_kpath .eq. 0) then
        deallocate(gloc)
        deallocate(projected_g)
        deallocate(delta_mat)
        deallocate(e_projected_mat)
        if (cal_mode .eq. 2) then ! calclating DOS after analytic continuation
          deallocate(tdos)
          deallocate(pdos)
        endif
      elseif (cal_mode .eq. 3) then ! calculating A(k,omega) after analytic continuation
        deallocate(g_path_trace)
        deallocate(g_path_orb)
        deallocate(g_path)
      endif

      end


! cal_h_qp never writes out wannier_den_matrix. It also pays no attention to cal_mode.
      subroutine cal_h_qp(numk,ndim_k,n_mpi_k,ndim_k_max,kvec)

      use comlowh_mod
      use pfit_mod

      implicit none
      include "mpif.h"

      integer, intent(in) :: numk, ndim_k(nproc), n_mpi_k(nproc),
     $  ndim_k_max
      double precision, intent(in) :: kvec(3,numk)

      integer ::
     $  ii,jj,i1,j1,ik,irk,iomega,ibnd,jbnd,iwancell,jwancell,kwancell,
     $  iorb,jorb,iatom,icell,jcell,iwan,jwan,igshift,
     $  iio,fail_ind,proc_num,indout,tempind(1)

      double precision ::pi,inv_numk,
     $  deltaij,diag(num_wancell),
     $  coeff1(4),
     $  coeff2(4),      
     $  zinv_m1_imp(n_ioac,n_iac_mat_i:n_iac_mat_f),        
     $  sig0(n_ioac,n_iac_mat_i:n_iac_mat_f),
     $  tempvec(num_wancell)
      complex*16 :: ai, phaseval_dcmplx(num_cell),
     $  h_k_wan_co(nqdiv),
     $  h_r_wan(nrpts),
     $  tempmat1(num_orb_max,num_wancell),
     $  tempmat2(num_orb_max,num_orb_max),
     $  tempmat3(num_wancell,num_wancell),
     $  zinv(num_wancell,num_wancell),
     $  zinv_imp(num_wancell,num_wancell),      
     $  zinv_half(num_wancell,num_wancell),
     $  z_half_imp(num_wancell,num_wancell),

     $  mat_loc_ackw(num_orb_max,num_orb_max),
     $  tempmat5(num_orb_max,num_orb_max),
     $  tempmat11(num_wancell, num_wancell),
     $  tempmat12(num_wancell, num_wancell),
     $  fac,
     $  h_qp_tmp(num_wancell,num_wancell)


      integer :: irequest,istatus(mpi_status_size)      
      logical :: mpiflag
      logical :: file_exist      
      
      double precision, allocatable :: tempdble(:,:),
     $  eigenvalues_fi_all(:,:),
     $  eigenvalues_qp_fi_all(:,:),
     $  eigenvalues_qp_orb(:,:),
     $  eigenvalues_orb(:,:)            
      complex*16, allocatable :: tempmat4(:,:),h_qp(:,:,:)


c$$$      write(179,*) 'here1'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     1. zone_folding information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ai=dcmplx(0.0d0, 1.0d0)
      pi=datan2(1.0d0, 1.0d0)*4.0d0

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     2. h_k in fine grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      allocate(h_k_wan_fi(num_wancell,num_wancell,ndim_k_max))
      h_k_wan_fi=0.0d0

      do iwan=1, num_wann
        do jwan=1, num_wann
          
          h_k_wan_co=0.0d0
          
          do irk=1, nqdiv
            do ibnd=1, num_bands                
              h_k_wan_co(irk)=h_k_wan_co(irk)+
     $          eigenvalues(ibnd,irk)*v_matrix(ibnd,jwan,irk)
     $          *dconjg(v_matrix(ibnd,iwan,irk))
            enddo
          enddo
          
!     fourier transform to r space
          
          h_r_wan=0.0d0          
          do i1=1, nrpts
            do j1=1, nqdiv
              h_r_wan(i1)=h_r_wan(i1)
     $          +h_k_wan_co(j1)*
     $          cdexp(-2.0*ai*pi*
     $          sum(irvec(:,i1)*kpt_latt(:,j1)))
     $          /dble(nqdiv)
            enddo
          enddo

          if (num_cell .eq. 1) then ! num_cell=1
            do ik=1, ndim_k(me1)
              irk=n_mpi_k(me1)+ik
              
              do j1=1, nrpts
                h_k_wan_fi(iwan,jwan,ik)
     $            =h_k_wan_fi(iwan,jwan,ik)
     $            +h_r_wan(j1)*
     $            cdexp(2.0*ai*pi*
     $            sum(
     $            kvec(:,irk)
     $            *irvec(:,j1))
     $            )/dble(ndegen(j1))
              enddo
            enddo
c$$$          else
c$$$            do ik=1, ndim_k(me1)
c$$$              irk=n_mpi_k(me1)+ik
c$$$              do icell=1, num_cell
c$$$                do jcell=1, num_cell
c$$$                  iwancell=iwan+(icell-1)*num_wann
c$$$                  jwancell=jwan+(jcell-1)*num_wann
c$$$                  do igshift=1,num_cell ! gvec_conv_shift_fi vector
c$$$                    do j1=1, nrpts
c$$$                      phaseval
c$$$     $                  =sum(irvec(:,j1)
c$$$     $                  *matmul(
c$$$     $                  matmul(transpose(real_lattice),
c$$$     $                  recip_lattice_conv),
c$$$     $                  (kvec(:,irk)+gvec_conv_shift(:,igshift,irk))
c$$$     $                  )
c$$$     $                  )
c$$$                      phasetau
c$$$     $                  =sum((lpoint(:,icell)-lpoint(:,jcell))
c$$$     $                  *matmul(recip_lattice_conv,
c$$$     $                  (kvec(:,irk)+gvec_conv_shift(:,igshift,irk))
c$$$     $                  )
c$$$     $                  )
c$$$                      
c$$$                      h_k_wan_fi(iwancell,jwancell,ik)
c$$$     $                  =h_k_wan_fi(iwancell,jwancell,ik)
c$$$     $                  +h_r_wan(j1)*
c$$$     $                  cdexp(ai*(phaseval+phasetau))
c$$$     $                  /dble(ndegen(j1))/dble(num_cell)
c$$$                    enddo
c$$$                  enddo      
c$$$                enddo
c$$$              enddo    
c$$$            enddo
          endif
        enddo
      enddo
c$$$      write(179,*) 'here2'      
      
      allocate(eigenvalues_fi(num_wancell,ndim_k_max)) !
      eigenvalues_fi=0.0d0
      allocate(eigenvalues_fi_all(num_wancell,numk))
      eigenvalues_fi_all=0.0d0
      allocate(eigenvalues_qp_fi_all(num_wancell,numk))
      eigenvalues_qp_fi_all=0.0d0
      allocate(eigenvalues_orb(num_wancell,numk))
      eigenvalues_orb=0.0d0
      allocate(eigenvalues_qp_orb(num_wancell,numk))
      eigenvalues_qp_orb=0.0d0            
      allocate(w_matrix_fi(num_wancell,num_wancell,ndim_k_max))
      w_matrix_fi=0.0d0
      
      w_matrix_fi=h_k_wan_fi

! this is a file read
      eigenvalues_fi_all=0.0d0
      if (maswrk) then
        open(unit=10, file='ef.dat', status='old')
        read(10,*) ef
        close(10)
      endif

      
       mpi_bcast(ef,1,mpi_double_precision,0,
     $  MPI_COMM_WORLD,mpierr)

      do ik=1, ndim_k(me1)
        irk=n_mpi_k(me1)+ik
        call hermitianeigen_cmplxdouble
     $    (num_wancell,eigenvalues_fi(1,ik),w_matrix_fi(1,1,ik))
        eigenvalues_fi_all(:,irk)=eigenvalues_fi(:,ik)
        do iwancell=1, num_wancell
          tempvec=cdabs(w_matrix_fi(iwancell,:,ik))
          tempind=maxloc(tempvec)
          eigenvalues_orb(iwancell,irk)=eigenvalues_fi(tempind(1),ik)
        enddo
        
        
      enddo

      call mpi_allreduce_dble(eigenvalues_fi_all,
     $  numk*num_wancell,mpi_sum,mpi_comm_world)
      call mpi_allreduce_dble(eigenvalues_orb,
     $  numk*num_wancell,mpi_sum,mpi_comm_world)      
c$$$  write(179,*) 'here3'
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     renormalize projector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      allocate(proj
     $  (num_orb_max,num_wancell,num_atom,num_cell,ndim_k_max))
      allocate(proj_renorm
     $  (num_orb_max,num_wancell,num_atom,num_cell,ndim_k_max))

      
      proj=0.0d0
      proj_renorm=0.0d0      

! dmft_projector changes the indexing of the input w_matrix_fi, which is num_wann x num_wann x ndim_k_max.  The output matrix proj is num_orb_max x num_wann x num_atom x num_cell=1 x ndim_k_max.  The first index of umat is transformed into the first and third indices of the output matrix.
! dmft_projector also applies a projector P_c on the left side of w_matrix_fi; the input matrix w_matrix_fi varies over all wanniers, and the output matrix proj contains only correlated wanniers.
! numk and kvec are not used.
      call dmft_projector
     $  (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,
     $  w_matrix_fi,proj)

      ! is_kpath means we are calculating qp spectra
      if (is_kpath .eq. 0) then

! cal_hmat_rhomat_from_proj calculates and prints out a hamiltonian
!   and a density matrix. It does not return anything to its caller.
        call cal_hmat_rhomat_from_proj
     $    (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,
     $    eigenvalues_fi,proj,
     $    ', bare projector, fine grid')

      endif

c$$$      write(179,*) 'here4'
        
      allocate(nbnd_in_win(numk))
      allocate(bnd_in_win(num_wancell,numk))      
      
! find_bnd_in_win creates bnd_in_win which says which bands are inside the energy window [ewin_min,ewin_max] for at least one k-point, and nbnd_in_win which is the count of bands inside the window.  Even though bnd_in_win and nbnd_in_win are indexed over k-points, they have the same value for all k-points.
      call find_bnd_in_win(eigenvalues_fi_all,numk,num_wancell,
     $  ewin_min,ewin_max,den_mode,is_spinorbit,nbnd_in_win,bnd_in_win)
c$$$      if (maswrk) then
c$$$        write(iun,*) 'num_band information in a find grid'
c$$$        write(iun,'(a, 10000i5)') 'bands in the window',
c$$$     $    maxval(nbnd_in_win), nbnd_in_win
c$$$        do irk=1, numk
c$$$          do ibnd=1, nbnd_in_win(irk)
c$$$            jbnd=bnd_in_win(ibnd,irk)
c$$$            write(iun,'(2i5, 2f12.6)')
c$$$     $        ibnd, jbnd,
c$$$     $        eigenvalues_fi_all(jbnd,irk)
c$$$          enddo
c$$$        enddo
c$$$      endif
      
! projector_renormalization_allk returns proj_out =  ( ( \sqrt{proj_in * P_E *  proj_in^\dagger})^{-1/2} ) * proj_in * P_E.  proj_out has the property that proj_out proj_out^\dagger = 1.
! P_E is a projector which keeps only the sector defined by ind_bnd_k, which says which bands are inside a an energy window.
! input: proj_in(num_orb_max,num_wann,num_atom,num_cell,ndim_kk_co_max). Structurally this means proj_in(orb+atom,wann,k)
! output: proj_out(num_orb_max,num_wann,num_atom,num_cell,ndim_kk_co_max) Structurally this means proj_out(orb+atom,wann,k)
      call projector_renormalization_allk
     $  (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,    
     $  eigenvalues_fi,proj,nbnd_in_win,bnd_in_win,proj_renorm)

      ! is_kpath means we are calculating qp spectra
      if (is_kpath .eq. 0) then

! cal_hmat_rhomat_from_proj calculates and prints out a hamiltonian
!   and a density matrix. It does not return anything to its caller.
        call cal_hmat_rhomat_from_proj
     $    (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,
     $    eigenvalues_fi,proj_renorm,
     $    ', renormalized projector, fine grid')

! projector_compare does nothing unless (is_kpath .eq. 1), which is not true here.
! projector_compare calculates proj_orig * proj_ren^\dagger and prints the result.
        call projector_compare
     $    (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,          
     $    proj,proj_renorm,'fine_grid')

      endif
c$$$      write(179,*) 'here5'
      

      inquire(file="trans_basis.dat", exist=is_trans_basis)
! In summary, trans_basis is equal to the identity and has no effect on the calculation.

      
      if (is_trans_basis) then      
        allocate(trans_basis
     $    (num_orb_max,num_orb_max,n_iac_mat_i:n_iac_mat_f))
        trans_basis=0.0d0

        if (maswrk) then
          open(unit=305, file='trans_basis.dat')
          allocate(tempdble(2,num_orb_max))
          do ii=n_iac_mat_i, n_iac_mat_f
            tempdble=0.0d0            
            if (ii .eq. 0) cycle
            do jorb=1, norb_iac(ii)
              read(305,*) ((tempdble(jj,iorb), jj=1, 2), 
     $          iorb=1,norb_iac(ii))
              do iorb=1, norb_iac(ii)
                trans_basis(jorb,iorb,ii)
     $            =dcmplx(tempdble(1,iorb), tempdble(2,iorb))
              enddo
            enddo
c$$$  test if it is unitary
            call check_if_unitary(norb_iac(ii),num_orb_max,
     $        trans_basis(1,1,ii),fail_ind)
            if (fail_ind .ne. 0) then
              write(iun,*) 'trans_basis matrix is not unitary'
              call ending                
            endif
          enddo
          deallocate(tempdble)
          close(305)
        endif
        call mpi_bcast(trans_basis,n_iac_mat*num_orb_max**2,
     $    mpi_double_complex,0,MPI_COMM_WORLD,mpierr)
      endif


      if (is_trans_basis) then
! In summary, trans_basis is equal to the identity and has no effect on the calculation.

      
        do icell=1, num_cell
          do iatom=1, num_atom
            do ik=1, ndim_k(me1)
              irk=n_mpi_k(me1)+ik                          
              ii=iac(iatom,icell)
              tempmat1=proj_renorm(:,:,iatom,icell,ik)
              proj_renorm(:,:,iatom,icell,ik)=0.0d0              
              call zgemm('n','n',num_orb(iatom),num_wancell,
     $          num_orb(iatom),(1.d0,0.d0),
     $          trans_basis(1,1,ii),num_orb_max,
     $          tempmat1,num_orb_max,
     $          (0.d0,0.d0),
     $          proj_renorm(1,1,iatom,icell,ik),num_orb_max)
            enddo
          enddo
        enddo
      endif

c$$$      write(179,*) 'here6'
      

c$$$  self_energy

! Here is is the origin of sig=sig.dat:
! comdmft.py's run_impurity_solver runs CTQMC and EVALSIM, and then:
! reads in ctqmc's output from params.obs.json and updates convergence.log
! "self-energy" from params.obs.json is saved in sig_bare.dat
! "self-energy" is smoothed using Gaussian broadening and stored in sigma.
!    sigma is saved in  sig_smth.dat.
! If any element of imag(sigma) is positive, then sig_causality is set
!       to 0=False; otherwise it is 1=True.
! If any element of imag(sigma) is positive, then sigma_to_delta = sigma_old [read in from sig.dat].
!    If it is true then sigma_to_delta is a mix of sigma with sigma_old [ read in from sig.dat].
! sigma_to_delta is saved in sig.dat.
! In summary, sig.dat comes from params.obs.json, then gaussian smoothing, then a check for negative imaginary parts, then mixing.
! ------------
! On the initial self-consistent cycle of dft+dmft and lqsgw+dmft, things are different:
! comdmft.py's generate_initial_self_energy creates sig.dat.
! If ('initial_self_energy' in control,  comdmft.py copies that self_energy to sig.dat.
!   Possibly it also copies information from initial_impurity_dir.
! Otherwise, comdmft.py copies data from dc.dat to sig.dat.
! ----------------
! If doing lqsgw+dmft, then before doing the self-consistent cycle there is a comlowh+comdc run.  To prepare for this dc.dat is all zeros, and sig.dat contains all zero's and omega's.  This is set up by comdmft.py's prepare_seed_dc_sig_and_wannier_dat.
      zinv_m1_imp=0.0d0
      sig0=0.0d0
      do ii=n_iac_mat_i, n_iac_mat_f
        if ((is_magnetic .eq. 1) .and. (ii .eq. 0)) cycle
        do iio=1, maxval(ioac(:,:,ii))
          call pfit(frequency_dble(1:6),dble(sig(iio,ii,1:6)),
     $      (/1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0/),
     $      (/0,1,2,3/), coeff1)          
          call pfit(frequency_dble(1:6),dimag(sig(iio,ii,1:6)),
     $      (/1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0/),
     $      (/0,1,2,3/), coeff2)
          sig0(iio,ii)=coeff1(1)          
          zinv_m1_imp(iio,ii)=-coeff2(2)
          if (maswrk) then
            write(iun,'(a, 2i5,6(2x,f12.6))') 'sig fit real part',
     $        ii,iio,coeff1
            write(iun,'(a, 2i5,6(2x,f12.6))') 'sig fit imag part',
     $        ii,iio,coeff2            
          endif
        enddo
      enddo
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     5. projected green's function at each frequency
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(h_qp(num_wancell,num_wancell,ndim_k_max))
      h_qp=0.0d0
      
      allocate(h_nonint(num_wancell,num_wancell,ndim_k_max))
      h_nonint=0.0d0

      do ik=1, ndim_k(me1)
        irk=n_mpi_k(me1)+ik

        if (zinv_m1_exist) then

          zinv=0.0d0
          do iwancell=1, num_wancell
            zinv(iwancell,iwancell)=1.0d0
          enddo
          
          do icell=1, num_cell ! num_cell=1
            do iatom=1, num_atom
              ii=iac(iatom,icell)
              mat_loc_ackw=0.0d0
              do iorb=1, num_orb(iatom)
                do jorb=1, num_orb(iatom)
                  iio=ioac(iorb,jorb,ii)
                  if (iio .ne. 0) then
                    mat_loc_ackw(iorb,jorb)
     $                =mat_loc_ackw(iorb,jorb)
     $                +zinv_m1(iio,ii)
                  endif
                enddo
              enddo
              
              tempmat1=0.0d0
              call zgemm('n','n',num_orb(iatom),num_wancell,
     $          num_orb(iatom),(1.d0,0.d0),
     $          mat_loc_ackw,num_orb_max,
     $          proj_renorm(1,1,iatom,icell,ik),num_orb_max,
     $          (0.d0,0.d0),
     $          tempmat1,num_orb_max)
              
              call zgemm('c','n',num_wancell,num_wancell,
     $          num_orb(iatom),(1.d0,0.d0),
     $          proj_renorm(1,1,iatom,icell,ik),num_orb_max,
     $          tempmat1,num_orb_max,(1.d0,0.d0),
     $          zinv,num_wancell)
            enddo
          enddo
          
          call hermitianeigen_cmplxdouble
     $      (num_wancell,diag,zinv)
          
          zinv_half=0.0d0
          do iwancell=1, num_wancell
            do jwancell=1, num_wancell
              do kwancell=1, num_wancell                              
                zinv_half(iwancell,jwancell)
     $            =zinv_half(iwancell,jwancell)
     $            +zinv(iwancell,kwancell)*dsqrt(diag(kwancell))
     $            *dconjg(zinv(jwancell,kwancell))
              enddo
            enddo
          enddo
          
          do iwancell=1, num_wancell
            do jwancell=1, num_wancell
              do kwancell=1, num_wancell                              
                h_nonint(iwancell,jwancell,ik)
     $            =h_nonint(iwancell,jwancell,ik)
     $            +zinv_half(iwancell,kwancell)
     $            *(eigenvalues_fi(kwancell,ik))
     $            *zinv_half(kwancell,jwancell)
              enddo
            enddo
          enddo

        else ! not (zinv_m1_exist)

          do iwancell=1, num_wancell
            h_nonint(iwancell,iwancell,ik)
     $        =eigenvalues_fi(iwancell,ik)
          enddo

        endif ! not (zinv_m1_exist)

        do iwancell=1,num_wancell
          h_nonint(iwancell,iwancell,ik)
     $      =h_nonint(iwancell,iwancell,ik)-ef
        enddo

      enddo ! ik=1, ndim_k(me1)


        
c$$$      write(179,*) 'here7'
        
      
! This is the origin of dc.dat: If doing qsgw+dmft, then run_dc runs ComDC, and then takes data (dc[0,jj,kk]) from ComDC's sig_mat.dat and puts it dc_mat.dat.  Alternatively, if doing dft+dmft, then run_dc/cal_nominal_dc creates dc_mat.dat.  Alternatively, if ('dc_mat_to_read' in control), then run_dc only copies an old version of dc_mat.dat from somewhere and uses that. After ComLowH is run, delta_postprocessing/cal_dc_diagonal takes takes from dc_mat.dat and puts it in dc.dat.
! when doing lqsgw+dmft, during the first call to comlowh + comdc (before the loop starts), dc.dat is all zeros.  This is set up by comdmft.py's prepare_seed_dc_sig_and_wannier_dat.
! dc.dat is used both by ComLowH and by CTQMC.
      eigenvalues_qp_fi_all=0.0d0     

      do ik=1, ndim_k(me1)
        irk=n_mpi_k(me1)+ik        
        
!     embeding self-energy term:E(Sigma-DC)
        
        do icell=1, num_cell ! num_cell=1
          do iatom=1, num_atom

            ii=iac(iatom,icell)
            mat_loc_ackw=0.0d0

            do iorb=1, num_orb(iatom)
              do jorb=1, num_orb(iatom)

                iio=ioac(iorb,jorb,ii)
                if (iio .ne. 0) then
                  mat_loc_ackw(iorb,jorb)
     $              =mat_loc_ackw(iorb,jorb)
     $              +sig0(iio,ii)
     $              -dc(iio,ii)
                endif
                
              enddo ! jorb=1, num_orb(iatom)
            enddo ! iorb=1, num_orb(iatom)
            
            tempmat1=0.0d0
            call zgemm('n','n',num_orb(iatom),num_wancell,
     $        num_orb(iatom),(1.d0,0.d0),
     $        mat_loc_ackw,num_orb_max,
     $        proj_renorm(1,1,iatom,icell,ik),num_orb_max,
     $        (0.d0,0.d0),
     $        tempmat1,num_orb_max)
            
            call zgemm('c','n',num_wancell,num_wancell,
     $        num_orb(iatom),(1.d0,0.d0),
     $        proj_renorm(1,1,iatom,icell,ik),num_orb_max,
     $        tempmat1,num_orb_max,(1.d0,0.d0),
     $        h_nonint(1,1,ik),num_wancell)

          enddo ! iatom=1, num_atom
        enddo  ! icell=1, num_cell ! num_cell=1


        zinv_imp=0.0d0
        do iwancell=1, num_wancell
          zinv_imp(iwancell,iwancell)=1.0d0
        enddo
        
        do icell=1, num_cell ! num_cell=1
          do iatom=1, num_atom

            ii=iac(iatom,icell)
            mat_loc_ackw=0.0d0
            do iorb=1, num_orb(iatom)
              do jorb=1, num_orb(iatom)
                iio=ioac(iorb,jorb,ii)
                if (iio .ne. 0) then
                  mat_loc_ackw(iorb,jorb)
     $              =mat_loc_ackw(iorb,jorb)
     $              +zinv_m1_imp(iio,ii)
                endif
              enddo
            enddo
            
            tempmat1=0.0d0
            call zgemm('n','n',num_orb(iatom),num_wancell,
     $        num_orb(iatom),(1.d0,0.d0),
     $        mat_loc_ackw,num_orb_max,
     $        proj_renorm(1,1,iatom,icell,ik),num_orb_max,
     $        (0.d0,0.d0),
     $        tempmat1,num_orb_max)
            
            call zgemm('c','n',num_wancell,num_wancell,
     $        num_orb(iatom),(1.d0,0.d0),
     $        proj_renorm(1,1,iatom,icell,ik),num_orb_max,
     $        tempmat1,num_orb_max,(1.d0,0.d0),
     $        zinv_imp,num_wancell)

          enddo ! iatom=1, num_atom
        enddo ! icell=1, num_cell ! num_cell=1
        
        call hermitianeigen_cmplxdouble
     $    (num_wancell,diag,zinv_imp)
        z_half_imp=0.0d0
        do iwancell=1, num_wancell
          do jwancell=1, num_wancell
            do kwancell=1, num_wancell                              
              z_half_imp(iwancell,jwancell)
     $          =z_half_imp(iwancell,jwancell)
     $          +zinv_imp(iwancell,kwancell)
     $          *1.0d0/dsqrt(diag(kwancell))
     $          *dconjg(zinv_imp(jwancell,kwancell))
            enddo
          enddo
        enddo
        
        
        call zgemm('n','n',num_wancell,num_wancell,num_wancell,
     $    (1.0d0,0.0d0), z_half_imp(1,1),num_wancell,
     $    h_nonint(1,1,ik),num_wancell,          
     $    (0.d0,0.d0),tempmat11,num_wancell)
        
        call zgemm('n','n',num_wancell,num_wancell,num_wancell,
     $    (1.d0,0.d0),tempmat11,num_wancell,
     $    z_half_imp(1,1),num_wancell,
     $    (0.d0,0.d0),tempmat12,num_wancell)
        
        call zgemm('n','n',num_wancell,num_wancell,num_wancell,
     $    (1.0d0,0.0d0), w_matrix_fi(1,1,ik),num_wancell,
     $    tempmat12,num_wancell,          
     $    (0.d0,0.d0),tempmat11,num_wancell)
        
        call zgemm('n','c',num_wancell,num_wancell,num_wancell,
     $    (1.d0,0.d0),tempmat11,num_wancell,
     $    w_matrix_fi(1,1,ik),num_wancell,
     $    (0.d0,0.d0),h_qp(1,1,ik),num_wancell)
        
c$$$        do iwancell=1, num_wancell
c$$$          h_qp(iwancell,iwancell,ik)
c$$$     $      =h_qp(iwancell,iwancell,ik)-ef
c$$$        enddo


        tempmat12=h_qp(:,:,ik)
        call hermitianeigen_cmplxdouble
     $    (num_wancell,eigenvalues_qp_fi_all(1,irk),
     $    tempmat12)
c$$$  write(179,*) 'here8',irk,ik
       do iwancell=1, num_wancell
          tempvec=cdabs(tempmat12(iwancell,:))
          tempind=maxloc(tempvec)
          eigenvalues_qp_orb(iwancell,irk)
     $      =eigenvalues_qp_fi_all(tempind(1),irk)
        enddo
        
      enddo ! ik=1, ndim_k(me1)

c$$$      write(179,*) 'here8'
        

      
      call mpi_allreduce_dble(eigenvalues_qp_fi_all,
     $  numk*num_wancell,mpi_sum,mpi_comm_world)
c$$$  write(179,*) 'here9'
      call mpi_allreduce_dble(eigenvalues_qp_orb,
     $  numk*num_wancell,mpi_sum,mpi_comm_world)
      
      
      if (maswrk) then
        open(unit = 990, file = 'wannier_band_qp_interpolated.dat')
        do ibnd=1,num_wancell
          do irk=1, numk
            write(990,'(2i7, f20.12)')
     $        irk, ibnd,eigenvalues_qp_fi_all(ibnd,irk)
          enddo
        enddo
        close(990)

        open(unit = 990, file = 'wannier_band_non_interpolated.dat')
        do ibnd=1,num_wancell
          do irk=1, numk
            write(990,'(2i7, f20.12)')
     $        irk, ibnd,eigenvalues_fi_all(ibnd,irk)
          enddo
        enddo
        close(990)

        open(unit = 990, file = 'wannier_band_qp_interpolated_orb.dat')
        do ibnd=1,num_wancell
          do irk=1, numk
            write(990,'(2i7, f20.12)')
     $        irk, ibnd,eigenvalues_qp_orb(ibnd,irk)
          enddo
        enddo
        close(990)

        open(unit = 990, file = 'wannier_band_non_interpolated_orb.dat')
        do ibnd=1,num_wancell
          do irk=1, numk
            write(990,'(2i7, f20.12)')
     $        irk, ibnd,eigenvalues_orb(ibnd,irk)
          enddo
        enddo
        close(990)        
      endif        


c$$$      !!!!!!!!!!!!!!!!!!! here
c$$$      if (maswrk) then      
c$$$        open(unit=990, file='h_qp_matrix.dat')
c$$$      endif
c$$$      
c$$$      if (.not. maswrk) then
c$$$        do ik=1, ndim_k(me1)
c$$$          irk=n_mpi_k(me1)+ik
c$$$          call mpi_isend(h_qp(1,1,ik),
c$$$     $      num_wancell**2,mpi_double_complex,0,irk,
c$$$     $      mpi_comm_world,irequest,mpierr)
c$$$          call mpi_wait(irequest,istatus,mpierr)
c$$$          call mpi_test(irequest,mpiflag,istatus,mpierr)
c$$$          if (.not.mpiflag)
c$$$     $      call mpi_abort(mpi_comm_world,mpierr)
c$$$        enddo
c$$$      else
c$$$        do irk=1, numk
c$$$          call find_proc_ind(ndim_k, n_mpi_k, nproc,irk,
c$$$     $      proc_num,indout)
c$$$          if (proc_num .eq. 0)then
c$$$            write(990,'(2000(2f20.12,2x))')
c$$$     $        ((h_qp(iwancell,jwancell,indout),
c$$$     $        iwancell=1,num_wancell),
c$$$     $        jwancell=1,num_wancell)
c$$$          else
c$$$            h_qp_tmp=0.0d0
c$$$            call mpi_irecv(h_qp_tmp,
c$$$     $        num_wancell**2,mpi_double_complex,proc_num,
c$$$     $        irk,mpi_comm_world,irequest,mpierr)
c$$$            call mpi_wait(irequest,istatus,mpierr)
c$$$            call mpi_test(irequest,mpiflag,istatus,mpierr)
c$$$            if (.not.mpiflag)
c$$$     $        call mpi_abort(mpi_comm_world,mpierr)
c$$$            write(990,'(2000(2f20.12,2x))')
c$$$     $        ((h_qp_tmp(iwancell,jwancell),
c$$$     $        iwancell=1,num_wancell),
c$$$     $        jwancell=1,num_wancell)              
c$$$          endif
c$$$        enddo
c$$$      endif
c$$$      if (maswrk) then
c$$$        close(990)
c$$$      endif
c$$$      
c$$$      !!!!!!!!!!!!!!!here
c$$$      


      deallocate(nbnd_in_win) !
      deallocate(bnd_in_win) !
      
      deallocate(h_k_wan_fi) !
      deallocate(eigenvalues_fi) !
      deallocate(eigenvalues_fi_all) !
      deallocate(eigenvalues_qp_fi_all) !      
      deallocate(w_matrix_fi) !
      deallocate(h_nonint) !
      deallocate(proj)!
      deallocate(proj_renorm) !
      if (is_trans_basis) then            
        deallocate(trans_basis) !
      endif
      
      end
      


      subroutine size_shift_par(n,npr,ndim_pr,n_mpi_pr)

c$$$  parallization of processor
c$$$  among the "npr" processor, we allocate "n" jobs
c$$$  to illustrate if we have 5 processor and 8 jobs
c$$$  the first 3 processors work on 2 jobs and the last work on 1 job      
c$$$  ndim_pr: the number of jobs in each processor: ex) (2, 2, 2, 1, 1)
c$$$  n_mpi_pr: (the index of first jobs-1) in each processor ex) (0, 2, 4, 6, 7)

      implicit none
      integer, intent(in) :: n,npr
      integer, intent(out) :: ndim_pr(npr),n_mpi_pr(npr)
      integer :: ish,irest,ind,i
      ish=n/npr
      irest=n-ish*npr
      ndim_pr(1:irest)=ish+1
      ndim_pr(irest+1:npr)=ish
      ind=0
      do i=1,npr
        n_mpi_pr(i)=ind
        ind=ind+ndim_pr(i)
      enddo
      end

      SUBROUTINE mpi_allreduce_int(X,LENX,OP,comm)
      implicit none
      INCLUDE 'mpif.h'
      integer,intent(in) :: lenx,op,comm
      integer, intent(inout) :: x(lenx)
      integer :: mxbuff
      integer,allocatable :: buff(:)
      integer :: i,j,npass,length,locx,ierror

      mxbuff=8192
      allocate(buff(mxbuff))
      NPASS = (LENX-1)/MXBUFF + 1
      LENGTH = MXBUFF
      LOCX = 1
      DO I=1,NPASS
        IF(I.EQ.NPASS) LENGTH = LENX - MXBUFF*(NPASS-1)
        CALL MPI_ALLREDUCE(X(LOCX),BUFF,LENGTH,MPI_INTEGER,
     $    OP,comm,IERROR)
c$$$        if (ierror .ne. MPI_SUCCESS) then
c$$$          write(179,*) 'mpi_allreduce_int fail',ierror
c$$$          call flush(179)
c$$$        endif
c$$$  CALL ICOPY(LENGTH,BUFF,1,X(LOCX),1)
        do j=1,length
          x(locx+j-1)=buff(j)
        enddo
        LOCX = LOCX + MXBUFF
      enddo
      deallocate(buff)
      END

      SUBROUTINE mpi_allreduce_dble(X,LENX,OP,comm)
      implicit none
      INCLUDE 'mpif.h'
      integer,intent(in) :: lenx,op,comm
      double precision, intent(inout) :: x(lenx)
      integer :: mxbuff
      double precision, allocatable :: buff(:)
      integer :: i,j,npass,length,locx,ierror

      mxbuff=4096
      allocate(buff(mxbuff))
      NPASS = (LENX-1)/MXBUFF + 1
      LENGTH = MXBUFF
      LOCX = 1
      DO I=1,NPASS
        IF(I.EQ.NPASS) LENGTH = LENX - MXBUFF*(NPASS-1)
        CALL MPI_ALLREDUCE(X(LOCX),BUFF,LENGTH,MPI_double_precision,
     $    OP,comm,IERROR)
c$$$        if (ierror .ne. MPI_SUCCESS) then
c$$$          write(179,*) 'mpi_allreduce_dble fail',ierror
c$$$          call flush(179)
c$$$        endif
        CALL DCOPY(LENGTH,BUFF,1,X(LOCX),1)
c$$$  do j=1,length
c$$$  x(locx+j-1)=buff(j)
c$$$  enddo
        LOCX = LOCX + MXBUFF
      enddo
      deallocate(buff)
      END


      SUBROUTINE mpi_allreduce_dcmplx(X,LENX,OP,comm)
      implicit none
      INCLUDE 'mpif.h'
      integer,intent(in) :: lenx,op,comm
      complex*16, intent(inout) :: x(lenx)
      integer :: mxbuff
      complex*16,allocatable :: buff(:)
      integer :: i,j,npass,length,locx,ierror


      mxbuff=2048
      allocate(buff(mxbuff))
      NPASS = (LENX-1)/MXBUFF + 1
      LENGTH = MXBUFF
      LOCX = 1
      DO I=1,NPASS
        IF(I.EQ.NPASS) LENGTH = LENX - MXBUFF*(NPASS-1)
        CALL MPI_ALLREDUCE(X(LOCX),BUFF,LENGTH,MPI_double_complex,
     $    OP,comm,IERROR)
c$$$        if (ierror .ne. MPI_SUCCESS) then
c$$$          write(179,*) 'mpi_allreduce_dble fail',ierror
c$$$          call flush(179)
c$$$        endif
        CALL ZCOPY(LENGTH,BUFF,1,X(LOCX),1)
c$$$  do j=1,length
c$$$  x(locx+j-1)=buff(j)
c$$$  enddo
        LOCX = LOCX + MXBUFF
      enddo
      deallocate(buff)
      END


      subroutine dcmplx_matinv(mat, invmat, dim,dimmax)
      implicit none
      integer, intent(in) :: dim, dimmax
      complex*16, intent(in) :: mat(dimmax, dimmax)
      complex*16, intent(out) :: invmat(dimmax,dimmax)
      
      integer :: ipiv(dim),info
      complex*16 :: wrk(dim**2)

      invmat=mat
      call zgetrf(dim,dim,invmat,
     $  dimmax,ipiv,info)
c$$$      if (info .ne. 0) then
c$$$        write(179, *) 'error in zgetrf', info
c$$$      endif
      call zgetri(dim,invmat,dimmax,
     $  ipiv,wrk,dim**2,info)          
c$$$      if (info .ne. 0) then
c$$$        write(179, *) 'error in zgetri', info
c$$$      endif
      end


! This routine finds a Fermi energy ef (used in gmat) that solves the equation 0=cal_tos_analytical_tail_fun.
! cal_tos_analytical_tail_fun=Tr(nimp(b))-nelectron
! nimp(b) = -2 \beta^{-1} N_k^{-1} \sum_k \sum_{\omega} \exp(-\beta \omega) diag[ gmat(k,\omega) - \omega^{-1} - g2mat(k,\omega_{max}) / \omega^2 - g3mat(k,\omega_{max}) / \omega^3] - N_k^{-1} \sum_k (-1/2 + g2mat(k,\omega_{max}) \beta/4)
! g2mat(k,\omega_{max}) = (g(k,\omega_{max}) + g^\dagger(k,\omega_{max}))/2 * \omega_{max}^2
! g3mat(k,\omega_{max})  = (g(k,\omega_{max}) - g^\dagger(k,\omega_{max})  - 2/\omega_{max})/2 * \omega_{max}^3
! gmat(k,\omega) = g(k,\omega) = (ef + frequency(\omega) - h_nonint - \sum_atoms proj_renorm^\dagger (sig(\omega) - dc - I * gamma_selfenergy) proj_renorm)^{-1}
! ---------
! g_gw is set to zero by the caller
! if no spin-orbit, nelectron=dble(nint(nelectron_dble*2))/2.0d0  +delta_charge/2.0d0
! if spin-orbit, nelectron=dble(nint(nelectron_dble))+delta_charge
! nelectron_dble=N_q^{-1} \sum_{q,b,o, eigenvalues(b,k) < 0} |v_matrix(b,o,k)|^2 
      subroutine search_fermi_energy_brent
     $  (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,g_gw)
      
      use comlowh_mod

      implicit none

      integer, intent(in) :: numk, ndim_k(nproc), n_mpi_k(nproc),
     $  ndim_k_max,g_gw
      double precision, intent(in) :: kvec(3,numk)

      integer :: ii
      double precision :: tos, beta_dos, machep, t,aa,bb,tos_aa,tos_bb,
     $  n_low,n_high, tos_ef,mu_low,mu_high,tos_ef_gw,n_target,mu_target
      double precision, external :: cal_tos_gw_fun,
     $  cal_tos_analytical_tail_fun,find_fermi_level

      if (maswrk) write(iun, *) 'search fermi energy'

      aa=-5.0d0
      bb=5.0d0
      machep = epsilon ( machep )
      t = 1.0d-5

c$$$  check if we have different sign at aa and bb
      if (g_gw .eq. 1) then ! This never executes because g_gw=0.
        tos_aa=cal_tos_gw_fun
     $    (aa,numk,ndim_k,n_mpi_k,ndim_k_max,kvec,0.0d0)
        tos_bb=cal_tos_gw_fun
     $    (bb,numk,ndim_k,n_mpi_k,ndim_k_max,kvec,0.0d0)        
      else
        tos_aa=cal_tos_analytical_tail_fun
     $    (aa,numk,ndim_k,n_mpi_k,ndim_k_max,kvec,0.0d0)
        tos_bb=cal_tos_analytical_tail_fun
     $    (bb,numk,ndim_k,n_mpi_k,ndim_k_max,kvec,0.0d0)              
      endif

      if ((tos_aa .gt. 0.0d0) .or. (tos_bb .lt. 0.0d0)) then
        if (maswrk) then
          write(iun,*) "wrong initial setup in search fermi level"
        endif
        call ending
      endif
      
      n_target=0.0d0
      if (g_gw .eq. 1) then  ! This never executes because g_gw=0.
        mu_target=find_fermi_level(aa, bb, machep,t, cal_tos_gw_fun,
     $    numk,ndim_k,n_mpi_k,ndim_k_max,kvec,nproc,n_target)
      else
        mu_target=find_fermi_level(aa, bb, machep,t,
     $    cal_tos_analytical_tail_fun,
     $    numk,ndim_k,n_mpi_k,ndim_k_max,kvec,nproc,n_target)
      endif

      ef=mu_target
      
      if (g_gw .eq. 1) then  ! This never executes because g_gw=0.
        tos_ef=cal_tos_gw_fun
     $    (ef,numk,ndim_k,n_mpi_k,ndim_k_max,kvec,0.0d0)
      else
        tos_ef=cal_tos_analytical_tail_fun
     $    (ef,numk,ndim_k,n_mpi_k,ndim_k_max,kvec,0.0d0)
        tos_ef_gw=cal_tos_gw_fun
     $    (ef,numk,ndim_k,n_mpi_k,ndim_k_max,kvec,0.0d0)        
      endif

      if (maswrk) then
        if (g_gw .eq. 1) then  ! This never executes because g_gw=0.
          write(iun,*) 'ef', ef, tos_ef
        else
          write(iun,*) 'ef', ef, tos_ef, tos_ef_gw
        endif

      endif            


      end      


      double precision function cal_tos_gw_fun(ef_tmp,
     $  numk,ndim_k,n_mpi_k,ndim_k_max,kvec, nshift)
      
      use comlowh_mod

      implicit none
      include 'mpif.h'
      
      integer, intent(in) :: numk, ndim_k(nproc), n_mpi_k(nproc),
     $  ndim_k_max
      double precision, intent(in) :: kvec(3,numk), ef_tmp,nshift
      integer :: ik, irk, ibnd

      cal_tos_gw_fun=0.0d0

      do ik=1, ndim_k(me1)
        irk=n_mpi_k(me1)+ik
        do ibnd=1, num_wancell
          cal_tos_gw_fun=cal_tos_gw_fun+1.0d0
     $      /(dexp(beta*(eigenvalues_fi(ibnd,ik)-ef_tmp))+1.0d0)
     $      /dble(numk)
        enddo
      enddo

      call mpi_allreduce_dble(cal_tos_gw_fun,1,mpi_sum,mpi_comm_world)
      cal_tos_gw_fun=cal_tos_gw_fun-nelectron*num_cell-nshift
      end            
      
! cal_n_matrix is almost exactly the same as cal_n_lattice+cal_tos_analytical_tail_fun, except 
!   that the former is evaluated on the coarse k-point grid associated with ndiv and nqdiv,
!   while the latter is evaluated on the fine k-point grid associated with ndiv_fi and nqdiv_fi.
!   The former is used to produce the density matrix which is supplied to flapwmbpt, while
!   the latter is used to produce an electron count for adjusting the Fermi level.
! cal_tos_analytical_tail_fun=Tr(nimp(b))-nelectron-nshift
! nimp(b) = -2 \beta^{-1} N_k^{-1} \sum_k \sum_{\omega} \exp(-\beta \omega) diag[ gmat(k,\omega) - \omega^{-1} - g2mat(k,\omega_{max}) / \omega^2 - g3mat(k,\omega_{max}) / \omega^3] - N_k^{-1} \sum_k (-1/2 + g2mat(k,\omega_{max}) \beta/4)
! g2mat(k,\omega_{max}) = (g(k,\omega_{max}) + g^\dagger(k,\omega_{max}))/2 * \omega_{max}^2
! g3mat(k,\omega_{max})  = (g(k,\omega_{max}) - g^\dagger(k,\omega_{max})  - 2/\omega_{max})/2 * \omega_{max}^3
! gmat(k,\omega) = g(k,\omega) = (ef_tmp + frequency(\omega) - h_nonint - \sum_atoms proj_renorm^\dagger (sig(\omega) - dc - I * gamma_selfenergy) proj_renorm)^{-1}
      subroutine cal_n_lattice
     $  (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,ef_tmp,
     $  nimp)
      
      use comlowh_mod

      implicit none
      include 'mpif.h'
      
      integer, intent(in) :: numk, ndim_k(nproc), n_mpi_k(nproc),
     $  ndim_k_max
      double precision, intent(in) :: kvec(3,numk),ef_tmp
      double precision, intent(out) :: nimp(num_wancell)
c$$$  $  g1_corr(num_wancell),g2_corr(num_wancell),
c$$$  $  g3_corr(num_wancell)
      integer :: ik, irk, ibnd, iomega      
      complex*16 :: tempval,
     $  g2mat(num_wancell,num_wancell),
     $  g3mat(num_wancell,num_wancell),
     $  gmat(num_wancell,num_wancell),
     $  g_lat_inv_temp(num_wancell,num_wancell),
     $  one_matrix_f(num_wancell,num_wancell)
      
      nimp=0.0d0
c$$$  g1_corr=0.0d0
c$$$  g2_corr=0.0d0
c$$$  g3_corr=0.0d0            
      one_matrix_f=0.0d0
      do ibnd=1, num_wancell
        one_matrix_f(ibnd,ibnd)=1.0d0
      enddo
      
      do ik=1, ndim_k(me1)
        irk=n_mpi_k(me1)+ik

        gmat=0.0d0        
        g_lat_inv_temp=g_lat_inv(:,:,ik,n_omega)        
        
! g_lat_inv = frequency(\omega_{max}) - h_nonint - \sum_atoms proj_renorm^\dagger (sig(\omega_{max}) - dc - I * gamma_selfenergy) proj_renorm
! g_lat_inv_temp = ef_tmp + frequency(\omega_{max}) - h_nonint - \sum_atoms proj_renorm^\dagger (sig(\omega_{max}) - dc - I * gamma_selfenergy) proj_renorm
        do ibnd=1, num_wancell
          g_lat_inv_temp(ibnd,ibnd)
     $      =g_lat_inv_temp(ibnd,ibnd)+ef_tmp
        enddo

c$$$  write(179,*), iomega, ik, irk

! gmat(k,\omega_{max}) = g(k,\omega_{max}) = (ef_tmp + frequency(\omega_{max}) - h_nonint - \sum_atoms proj_renorm^\dagger (sig(\omega_{max}) - dc - I * gamma_selfenergy) proj_renorm)^{-1}
        call dcmplx_matinv(g_lat_inv_temp, gmat,
     $    num_wancell,num_wancell)
        
! g2mat(k,\omega_{max}) = (g(k,\omega_{max}) + g^\dagger(k,\omega_{max}))/2 * \omega_{max}^2
! g3mat(k,\omega_{max})  = (g(k,\omega_{max}) - g^\dagger(k,\omega_{max})  - 2/\omega_{max})/2 * \omega_{max}^3
        g2mat=(gmat+transpose(dconjg(gmat)))
     $    /2.0d0*frequency(n_omega)**2
        g3mat=(gmat-transpose(dconjg(gmat))
     $    -one_matrix_f*2.0/frequency(n_omega))
     $    /2.0d0*frequency(n_omega)**3                  
        


        do iomega=1, n_omega
          
          gmat=0.0d0        

          g_lat_inv_temp=g_lat_inv(:,:,ik,iomega)        
          
          do ibnd=1, num_wancell
            g_lat_inv_temp(ibnd,ibnd)
     $        =g_lat_inv_temp(ibnd,ibnd)+ef_tmp
          enddo

c$$$  write(179,*), iomega, ik, irk

          call dcmplx_matinv(g_lat_inv_temp, gmat,
     $      num_wancell,num_wancell)
! gmat(k,\omega) = g(k,\omega) = (ef_tmp + frequency(\omega) - h_nonint - \sum_atoms proj_renorm^\dagger (sig(\omega) - dc - I * gamma_selfenergy) proj_renorm)^{-1}

          do ibnd=1,num_wancell
            nimp(ibnd)
     $        =nimp(ibnd)
     $        -2.0d0/beta*dble(
     $        cdexp(-beta*frequency(iomega))
     $        *(
     $        gmat(ibnd,ibnd)
     $        -1.0d0/frequency(iomega)
     $        -g2mat(ibnd,ibnd)/(frequency(iomega))**2
     $        -g3mat(ibnd,ibnd)/(frequency(iomega))**3                  
     $        )
     $        )/dble(numk)
c$$$  g1_corr(ibnd)
c$$$  $        =g1_corr(ibnd)
c$$$  $        -2.0d0/beta*dble(
c$$$  $        cdexp(-beta*frequency(iomega))
c$$$  $        *(
c$$$  $        -1.0d0/frequency(iomega)
c$$$  $        )
c$$$  $        )/dble(numk)
c$$$  g2_corr(ibnd)
c$$$  $        =g2_corr(ibnd)
c$$$  $        -2.0d0/beta*dble(
c$$$  $        cdexp(-beta*frequency(iomega))
c$$$  $        *(
c$$$  $        -g2mat(ibnd,ibnd)/(frequency(iomega))**2
c$$$  $        )
c$$$  $        )/dble(numk)
c$$$  g3_corr(ibnd)
c$$$  $        =g3_corr(ibnd)
c$$$  $        -2.0d0/beta*dble(
c$$$  $        cdexp(-beta*frequency(iomega))
c$$$  $        *(
c$$$  $        -g3mat(ibnd,ibnd)/(frequency(iomega))**3                  
c$$$  $        )
c$$$  $        )/dble(numk)            
          enddo ! ibnd=1,num_wancell
        enddo ! iomega=1, n_omega

        do ibnd=1,num_wancell        
!     asymptotic interval
          nimp(ibnd)
     $      =nimp(ibnd)
     $      -1.0d0
     $      *(
     $      -1.0d0/2.0d0            
     $      +dble(g2mat(ibnd,ibnd))*beta/4.0d0
c$$$  $        *(beta/beta-1.0d0/2.0d0)
c$$$  $        -g3mat(ibnd,ibnd)*beta**2/4.0d0
c$$$  $        *((beta/beta)**2-(beta/beta))
     $      )/dble(numk)

c$$$  g1_corr(ibnd)
c$$$  $      =g1_corr(ibnd)
c$$$  $      -1.0d0
c$$$  $      *(
c$$$  $      -1.0d0/2.0d0            
c$$$  c$$$  $        *(beta/beta-1.0d0/2.0d0)
c$$$  c$$$  $        -g3mat(ibnd,ibnd)*beta**2/4.0d0
c$$$  c$$$  $        *((beta/beta)**2-(beta/beta))
c$$$  $      )/dble(numk)          
c$$$  g2_corr(ibnd)
c$$$  $      =g2_corr(ibnd)
c$$$  $      -1.0d0
c$$$  $      *(
c$$$  $      +g2mat(ibnd,ibnd)*beta/4.0d0
c$$$  c$$$  $        *(beta/beta-1.0d0/2.0d0)
c$$$  c$$$  $        -g3mat(ibnd,ibnd)*beta**2/4.0d0
c$$$  c$$$  $        *((beta/beta)**2-(beta/beta))
c$$$  $      )/dble(numk)
        enddo ! ibnd=1,num_wancell 

      enddo ! ik=1, ndim_k(me1)
      
      call mpi_allreduce_dble(nimp,num_wancell,mpi_sum,mpi_comm_world)

c$$$  call mpi_allreduce_dble(g1_corr,
c$$$  $  num_wancell,mpi_sum,mpi_comm_world)
c$$$  call mpi_allreduce_dble(g2_corr,
c$$$  $  num_wancell,mpi_sum,mpi_comm_world)
c$$$  call mpi_allreduce_dble(g3_corr,
c$$$  $  num_wancell,mpi_sum,mpi_comm_world)                  

      end


! cal_n_matrix is almost exactly the same as cal_n_lattice+cal_tos_analytical_tail_fun, except 
!   that the former is evaluated on the coarse k-point grid associated with ndiv and nqdiv,
!   while the latter is evaluated on the fine k-point grid associated with ndiv_fi and nqdiv_fi.
!   The former is used to produce the density matrix which is supplied to flapwmbpt, while
!   the latter is used to produce an electron count for adjusting the Fermi level.
! cal_tos_analytical_tail_fun=Tr(nimp(b))-nelectron-nshift
! nimp(b) = -2 \beta^{-1} N_k^{-1} \sum_k \sum_{\omega} \exp(-\beta \omega) diag[ gmat(k,\omega) - \omega^{-1} - g2mat(k,\omega_{max}) / \omega^2 - g3mat(k,\omega_{max}) / \omega^3] - N_k^{-1} \sum_k (-1/2 + g2mat(k,\omega_{max}) \beta/4)
! g2mat(k,\omega_{max}) = (g(k,\omega_{max}) + g^\dagger(k,\omega_{max}))/2 * \omega_{max}^2
! g3mat(k,\omega_{max})  = (g(k,\omega_{max}) - g^\dagger(k,\omega_{max})  - 2/\omega_{max})/2 * \omega_{max}^3
! gmat(k,\omega) = g(k,\omega) = (ef_tmp + frequency(\omega) - h_nonint - \sum_atoms proj_renorm^\dagger (sig(\omega) - dc - I * gamma_selfenergy) proj_renorm)^{-1}
      double precision function cal_tos_analytical_tail_fun(ef_tmp,
     $  numk,ndim_k,n_mpi_k,ndim_k_max,kvec,nshift)
      use comlowh_mod

      implicit none
      include 'mpif.h'
      
      integer, intent(in) :: numk, ndim_k(nproc), n_mpi_k(nproc),
     $  ndim_k_max
      double precision, intent(in) :: kvec(3,numk),ef_tmp, nshift
      integer :: ik, irk, ibnd, iomega,ii
      double precision :: nimp(num_wancell)
c$$$  , g1_corr(num_wancell),
c$$$  $  g2_corr(num_wancell), g3_corr(num_wancell),
c$$$  $  g4_corr(num_wancell), g5_corr(num_wancell)

      call cal_n_lattice
     $  (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,ef_tmp,
     $  nimp)      
c$$$      call cal_n_lattice_gw_eigen      
c$$$     $  (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,ef_tmp,
c$$$  $  nimp)
c$$$      call cal_n_lattice_e_inf
c$$$     $  (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,ef_tmp,
c$$$     $  nimp)      
      
      cal_tos_analytical_tail_fun=sum(nimp)-nelectron*num_cell-nshift
      end


      SUBROUTINE ENDING
      implicit none
      include "mpif.h"
      integer :: ierror

c$$$      write(179,*) 'call ending'
c$$$      call flush(179)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      CALL MPI_FINALIZE(IERROR)
      stop      
      END                                                               

! this is a file read
! reads in wannier.dat
      subroutine read_wann_dat
      use comlowh_mod

      implicit none
      integer :: ik,ibnd,iwan
      double precision :: nelectron_dble
      
      open(unit=10, file='wannier.dat', form='unformatted')
      read(10) real_lattice
      read(10) recip_lattice      
      read(10) num_bands
      read(10) num_wann
      read(10) ndiv
      nqdiv=ndiv(1)*ndiv(2)*ndiv(3)
      allocate(kpt_latt(3,nqdiv))
      
      read(10) kpt_latt
      allocate(include_bands(num_bands))
      read(10) include_bands
      allocate(v_matrix(num_bands,num_wann,nqdiv))      
      read(10) v_matrix
      allocate(eigenvalues(num_bands,nqdiv))
      read(10) eigenvalues
      read(10) nrpts
      allocate(irvec(3,nrpts))      
      read(10) irvec
      allocate(ndegen(nrpts))
      read(10) ndegen
      close(10)

      nelectron_dble=0.0
      do ik=1, nqdiv
        do ibnd=1, num_bands
          do iwan=1, num_wann
            if (eigenvalues(ibnd,ik) .lt. 0.0d0) then
              nelectron_dble=nelectron_dble+1.0d0/dble(nqdiv)*
     $          cdabs(v_matrix(ibnd,iwan,ik))**2
            endif
          enddo
        enddo
      enddo      

      if (is_spinorbit .eq. 1) then
        nelectron=dble(nint(nelectron_dble*num_cell))/num_cell
     $    +delta_charge ! delta_charge comes from the 'doping' entry in comdmft.ini
      else
        nelectron=dble(nint(nelectron_dble*num_cell*2))/2.0d0/num_cell
     $    +delta_charge/2.0d0
      endif
      write(iun,'(a, 3f12.6)')
     $  'nelectron in one cell', nelectron_dble, 
     $  delta_charge, nelectron      
      end
      
! this is a file read
      subroutine read_wann_inip
      
      use comlowh_mod

      implicit none
      include "mpif.h"      
      integer :: iwan,iorb,ii,jj,jorb, iatom, lval,mval1
      double precision :: tempvec1(3),tempvec2(3),mval2,ival
c$$$      double precision,allocatable :: tempdble1(:,:)

      if (is_spinorbit .eq. 1) then
        num_alimj=4
      else
        num_alimj=3
      endif
      allocate(nalimj(num_wann))
      do iwan=1, num_wann
        nalimj(iwan)=1
      enddo
      max_nalimj=maxval(nalimj)      
      open(unit=111, file=adjustl(trim(seed_name))//'.inip') ! initial projector
c$$$  if(l==0) 'S '
c$$$  if(l==1) ' Y   Z   X '
c$$$  if(l==2) ' XY  YZ  2ZZ-XX-YY  XZ  XX-YY '
c$$$  if(l==3) ' 3XXY  XYZ  4YZZ-XXY-YYY 2ZZZ-3XXZ-3YYZX ... '
      allocate(alimj(num_alimj,max_nalimj,num_wann))
      alimj=0.0d0
      allocate(alimj_coeff(max_nalimj,num_wann))
      alimj_coeff=0.0d0
      allocate(alimj_xaxis(3,max_nalimj,num_wann))
      alimj_coeff=0.0d0
      allocate(alimj_zaxis(3,max_nalimj,num_wann))
      alimj_coeff=0.0d0                        
c$$$  allocate(tempdble1(2,max_nalimj))
      do iwan=1, num_wann
c$$$        tempdble1=0.0d0
        if (is_spinorbit .eq. 0) then
          read(111, *) 
     $      iatom, lval, mval1,
c$$$     $      (tempdble1(ii,1), ii=1, 2),
     $      (alimj_xaxis(ii, 1,iwan), ii=1, 3),
     $      (alimj_zaxis(ii, 1,iwan), ii=1, 3)
          alimj(1,1,iwan)=iatom
          alimj(2,1,iwan)=lval
          alimj(3,1,iwan)=mval1+lval+1
        else
          read(111, *) 
     $      iatom, lval, ival,mval2,
c$$$     $      (tempdble1(ii,1), ii=1, 2),
     $      (alimj_xaxis(ii, 1,iwan), ii=1, 3),
     $      (alimj_zaxis(ii, 1,iwan), ii=1, 3)
          alimj(1,1,iwan)=iatom
          alimj(2,1,iwan)=lval
          alimj(3,1,iwan)=nint(ival*2.0d0)
          alimj(4,1,iwan)=nint(mval2*2.0d0)
        endif
        do iorb=1, nalimj(iwan)
          alimj_coeff(iorb,iwan)=dcmplx(1.0d0, 0.0d0)
c$$$          alimj_coeff(iorb,iwan)=dcmplx(tempdble1(1,iorb),
c$$$     $      tempdble1(2,iorb))          
        enddo
      enddo
      close(111)
c$$$      deallocate(tempdble1)
      end
      
      subroutine dcross(a,b,c)
**********************************************************
*     cross product (ax,ay,az)=(bx,by,bz)*(cx,cy,cz)      *
**********************************************************
      implicit none
      double precision, intent(in) :: b(3),c(3)
      double precision, intent(out) :: a(3)
      a=0.0d0
      a(1)=b(2)*c(3)-b(3)*c(2)
      a(2)=b(3)*c(1)-b(1)*c(3)
      a(3)=b(1)*c(2)-b(2)*c(1)
      end

      subroutine reciprocal_lattice_cal(a, b)
      implicit none
      double precision, intent(in) :: a(3,3)
      double precision, intent(out) :: b(3,3)
      double precision :: dbletemp(3),pi, crossdot

      pi=datan2(1.0d0,1.0d0)*4.0d0
      
      call dcross(dbletemp,a(:,2),a(:,3))
      
      crossdot=sum(dbletemp*a(:,1))

      call dcross(dbletemp, a(:,2),a(:,3))
      b(:,1)=dbletemp/crossdot*2.0d0*pi
      call dcross(dbletemp, a(:,3),a(:,1))
      b(:,2)=dbletemp/crossdot*2.0d0*pi
      call dcross(dbletemp, a(:,1),a(:,2))
      b(:,3)=dbletemp/crossdot*2.0d0*pi
      end



      subroutine vol_ratio_cal(a, b, vol_ratio)
      implicit none
      double precision, intent(in) :: a(3,3),b(3,3)
      double precision, intent(out) :: vol_ratio
      double precision :: dbletemp(3), vol1, vol2

      call dcross(dbletemp,a(:,2),a(:,3))
      vol1=sum(dbletemp*a(:,1))
      call dcross(dbletemp,b(:,2),b(:,3))
      vol2=sum(dbletemp*b(:,1))
      vol_ratio=vol2/vol1
      end


! dmft_projector changes the indexing of the input matrix umat, which is num_wann x num_wann x ndim_k_max.  The output matrix proj_out is num_orb_max x num_wann x num_atom x num_cell=1 x ndim_k_max.  The first index of umat is transformed into the first and third indices of the output matrix.
! dmft_projector also applies a projector P_c on the left side of umat; the input matrix umat varies over all wanniers, and the output matrix proj_out contains only correlated wanniers.
! The mapping is done with the aid of wan_cor_index, which maps (iorb,iatom,icell) to an index over wanniers.
! numk and kvec are not used.
! w_matrix_fi <=> umat, proj<=>proj_out
      subroutine dmft_projector
     $  (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,umat,proj_out)


      use comlowh_mod

      implicit none
      include "mpif.h"

      integer, intent(in) :: numk, ndim_k(nproc), n_mpi_k(nproc),
     $  ndim_k_max
      double precision, intent(in) :: kvec(3,numk)

      complex*16, intent(in) :: umat(num_wancell,num_wancell,ndim_k_max)
      
      complex*16, intent(out) :: proj_out(num_orb_max,num_wancell,
     $  num_atom,num_cell,ndim_k_max)
      
      integer :: 
     $  iatom, isort,icell,
     $  atom_count,iorb,jorb,ik,irk,iband, ii,jj, k1,k2,k3,kk
      

! wan_cor_index maps (iorb,iatom,icell) to an index over wanniers.
! for instance, for a single d-shell atom, wan_cor_index=wan_correlated might be [8,9,10,11,12], and its iatom and icell indices are constrained to 1.

!     calculate projector in a fine grid
      proj_out=0.0d0

      do icell=1, num_cell ! num_cell = 1
        do iatom=1, num_atom
          do iorb=1, num_orb(iatom)
! num_orb(iatom) is the number of correlated wanniers on this atom.
            do ik=1, ndim_k(me1)
              irk=n_mpi_k(me1)+ik
              do iband=1, num_wancell ! num_wancell = num_wann
                proj_out(iorb,iband,iatom,icell,ik)
     $            =proj_out(iorb,iband,iatom,icell,ik)
     $            +umat(
     $            wan_cor_index(iorb,iatom,icell),
     $            iband,ik)
              enddo
            enddo
          enddo
        enddo
      enddo
      
      end

! projector_renormalization_allk returns proj_out =  ( ( \sqrt{proj_in * P_E *  proj_in^\dagger})^{-1/2} ) * proj_in * P_E.  proj_out has the property that proj_out proj_out^\dagger = 1.
! P_E is a projector which keeps only the sector defined by ind_bnd_k, which says which bands are inside a an energy window.
! input: proj_in(num_orb_max,num_wann,num_atom,num_cell,ndim_kk_co_max). Structurally this means proj_in(orb+atom,wann,k)
! output: proj_out(num_orb_max,num_wann,num_atom,num_cell,ndim_kk_co_max) Structurally this means proj_out(orb+atom,wann,k)
      subroutine projector_renormalization_allk
     $  (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,eigval,proj_in,nbnd_k,
     $  ind_bnd_k,proj_out)
      use comlowh_mod
      
      implicit none

      integer, intent(in) :: numk, ndim_k(nproc), n_mpi_k(nproc),
     $  ndim_k_max,nbnd_k(numk),ind_bnd_k(num_wancell,numk)
      double precision, intent(in) ::
     $  kvec(3,numk),eigval(num_wancell,ndim_k_max)
      complex*16, intent(in) :: proj_in
     $  (num_orb_max,num_wancell,num_atom,num_cell,ndim_k_max)
      complex*16, intent(out) :: proj_out
     $  (num_orb_max,num_wancell,num_atom,num_cell,ndim_k_max)
      integer :: ik,irk,ibnd,iorb,jorb,iatom,ind_temp,icell,ii,jj
      complex*16 :: proj_k(num_proj, num_wancell),
     $  proj_k_temp(num_proj, num_wancell)

      proj_out=0.0d0

! This maps from proj_in(orb+atom,wann,k) to  proj_k(num_proj=sum(num_orb),wann) * P_E.  P_E is a projector which keeps only the sector defined by ind_bnd_k, which says which bands are inside a an energy window.
      do ik=1, ndim_k(me1)

        irk=n_mpi_k(me1)+ik
        proj_k=0.0d0
        proj_k_temp=0.0d0
        do ibnd=1, nbnd_k(irk)
          do icell=1, num_cell ! num_cell=1
            do iatom=1, num_atom
              do iorb=1, num_orb(iatom)
                ind_temp=iorb+sum(num_orb(1:(iatom-1)))
     $            +(icell-1)*sum(num_orb)                  
                proj_k(ind_temp,ind_bnd_k(ibnd,irk))
     $            =proj_k(ind_temp,ind_bnd_k(ibnd,irk))
     $            +proj_in(iorb,ind_bnd_k(ibnd,irk),iatom,icell,ik)                    
              enddo
            enddo
          enddo
        enddo

! projector_renormalization returns proj_k_temp =  ( ( \sqrt{proj_k proj_k^\dagger})^{-1/2} ) * proj_k.  proj_k_temp has the property that proj_k_temp proj_k_temp^\dagger = 1.
        call projector_renormalization_dagger
     $    (num_proj,num_wancell,proj_k,proj_k_temp)

! This maps from proj_k(num_proj=sum(num_orb),wann) to proj_out(orb+atom,wann,k).
        do icell=1, num_cell ! num_cell = 1
          do iatom=1, num_atom
            do iorb=1, num_orb(iatom)
              ind_temp=iorb+sum(num_orb(1:(iatom-1)))
     $          +(icell-1)*sum(num_orb)
              do ibnd=1, num_wancell
                proj_out(iorb,ibnd,iatom,icell,ik)
     $            =proj_k_temp(ind_temp,ibnd)
              enddo
            enddo
          enddo
        enddo
      enddo
      end

! projector_renormalization returns proj_out = proj_in * ( \sqrt{proj_in^\dagger proj_in})^{-1/2}.  proj_out has the property that proj_out^\dagger proj_out = 1.
      subroutine projector_renormalization
     $  (dim1,dim2,proj_in,proj_out)
!     symmetric orthonormalization
      
      implicit none

      integer, intent(in) :: dim1, dim2
      complex*16, intent(in) :: proj_in(dim1,dim2)
      complex*16, intent(out) :: proj_out(dim1,dim2)

      integer :: iorb
      double precision :: diag(dim2)
      complex*16 :: overlapmat(dim2,dim2),overlap_mhalf(dim2,dim2),
     $  tempmat(dim2,dim2)
      
      proj_out=0.0d0

      overlapmat=0.0d0
      overlap_mhalf=0.0d0
      diag=0.0d0
      tempmat=0.0d0
      
        ! overlapmat = proj_in^\dagger proj_in
      call zgemm('c','n',dim2,dim2,dim1,(1.d0,0.d0),proj_in,dim1,
     $  proj_in,dim1,(0.d0,0.d0),overlapmat,dim2)
      
        ! diag and overlapmat are the eigenvalues and eigenvectors of overlapmat.
        ! overlapmat diag overlapmat^\dagger = old overlapmat (or h.c. I'm not sure)
      call hermitianeigen_cmplxdouble
     $  (dim2,diag,overlapmat)

        ! tempmat = overlapmat/\sqrt{diag}
c     ------ Here we get T = O^(-1/2) -----------------------------------
      do iorb=1,dim2
        tempmat(:,iorb)=overlapmat(:,iorb)/sqrt(diag(iorb))
      enddo

        ! overlap_mhalf = overlapmat diag^{-1/2} overlapmat^\dagger
        ! overlap_mhalf = 1/\sqrt{proj_in^\dagger proj_in} (or h.c. I'm not sure)
      call zgemm('n','c',dim2,dim2,dim2,
     $  (1.d0,0.d0),tempmat,dim2,overlapmat,dim2,
     $  (0.d0,0.d0),overlap_mhalf,dim2)

        ! proj_out = proj_in / \sqrt{proj_in^\dagger proj_in}
c$$$  --------- Orthonormalized projector ----------------------
      call zgemm('n','n',dim1,dim2,dim2,(1.d0,0.d0),
     $  proj_in,dim1,overlap_mhalf,dim2,(0.d0,0.d0),
     $  proj_out,dim1)

      end

! projector_renormalization returns proj_out =  ( ( \sqrt{proj_in proj_in^\dagger})^{-1/2} ) * proj_in.  proj_out has the property that proj_out proj_out^\dagger = 1.
      subroutine projector_renormalization_dagger
     $  (dim1,dim2,proj_in,proj_out)
!     symmetric orthonormalization
      
      implicit none

      integer, intent(in) :: dim1, dim2
      complex*16, intent(in) :: proj_in(dim1,dim2)
      complex*16, intent(out) :: proj_out(dim1,dim2)

      complex*16 :: proj_in_temp(dim2,dim1),
     $  proj_out_temp(dim2,dim1)

      proj_in_temp=dconjg(transpose(proj_in))
      
      call projector_renormalization
     $  (dim2,dim1,proj_in_temp,proj_out_temp)

      proj_out=dconjg(transpose(proj_out_temp))
      end
      

! projector_compare does nothing unless (is_kpath .eq. 1)
! projector_compare calculates proj_orig * proj_ren^\dagger and prints the result.
      subroutine projector_compare
     $  (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,proj_orig,proj_ren,
     $  textmessage)
      use comlowh_mod

      implicit none
      include "mpif.h"
      integer, intent(in) :: numk, ndim_k(nproc), n_mpi_k(nproc),
     $  ndim_k_max
      double precision, intent(in) :: kvec(3,numk)
      complex*16,intent(in) ::
     $  proj_orig(num_orb_max,num_wancell,num_atom,num_cell,ndim_k_max),
     $  proj_ren(num_orb_max,num_wancell,num_atom,num_cell,ndim_k_max)
      character(len=*) :: textmessage

      integer :: ik,irk,iorb,ibnd,jorb,iprj,jprj,iatom,icell,
     $  jcell,jatom,
     $  indtmp,jndtmp,iband,ii
      complex*16 :: overlap_mat(num_proj, num_proj),
     $  tempmat1_k(num_orb_max,num_atom,num_cell,num_wancell),
     $  tempmat2_k(num_orb_max,num_atom,num_cell,num_wancell),
     $  rhomat(num_orb_max,num_orb_max,num_atom,num_cell),
     $  hmat(num_orb_max,num_orb_max,num_atom,num_cell)      
      
      if (is_kpath .eq. 1) return
      
      overlap_mat=0.0d0

      do ik=1, ndim_k(me1)
        irk=n_mpi_k(me1)+ik
        tempmat1_k=0.0d0
        tempmat2_k=0.0d0
        do icell=1,num_cell
          do iatom=1, num_atom
            do iorb=1, num_orb(iatom)
              do ibnd=1, num_wancell
                tempmat1_k(iorb,iatom,icell,ibnd)
     $            =proj_orig(iorb,ibnd,iatom,icell,ik)
                tempmat2_k(iorb,iatom,icell,ibnd)
     $            =proj_ren(iorb,ibnd,iatom,icell,ik)
              enddo
            enddo
          enddo
        enddo
        
        do icell=1, num_cell
          do iatom=1, num_atom
            do iorb=1, num_orb(iatom)
              indtmp=iorb+sum(num_orb(1:(iatom-1)))
     $          +(icell-1)*sum(num_orb)              
              do jcell=1, num_cell
                do jatom=1, num_atom
                  do jorb=1, num_orb(jatom)
                    jndtmp=jorb+sum(num_orb(1:(jatom-1)))
     $                +(jcell-1)*sum(num_orb)                                  
                    do ibnd=1, num_wancell
                      overlap_mat(indtmp,jndtmp)
     $                  =overlap_mat(indtmp,jndtmp)
     $                  +tempmat1_k(iorb,iatom,icell,ibnd)
     $                  *dconjg(tempmat2_k(jorb,jatom,jcell,ibnd))
     $                  /dble(numk)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo        
      enddo      
      
      call mpi_allreduce_dcmplx(overlap_mat,(num_proj)**2,
     $  mpi_sum,mpi_comm_world)
      
      if (maswrk) then
        write(iun,*)
     $    "overlap_mat of new projector and old projector, "
     $    //textmessage
        do iprj=1, num_proj
          write(iun,'(100(2f12.6,2x))')
     $      (overlap_mat(iprj,jprj), jprj=1, num_proj)
        enddo
      endif
      
      end

      subroutine cal_gvec_conv_shift(kvec, gvec)
      use comlowh_mod
      implicit none
      include "mpif.h"      
      double precision,intent(in) :: kvec(3)
      integer, intent(out) :: gvec(3,num_cell)
      integer :: ia,ib,ic,ja,jb,jc,tempcnt,ii,notzero,
     $  issame(3,0:(num_cell-1),0:(num_cell-1),0:(num_cell-1)),
     $  degennum(0:(num_cell-1),0:(num_cell-1),0:(num_cell-1))
      double precision ::pi,tempdble3(3),
     $  dblearray(3,0:(num_cell-1),0:(num_cell-1),0:(num_cell-1))

      pi=datan2(1.0d0,1.0d0)*4.0d0

      dblearray=0.0d0
      issame=0
      degennum=0
      do ia=0, num_cell-1
        do ib=0, num_cell-1
          do ic=0, num_cell-1
            tempdble3=matmul(transpose(real_lattice)/2.0d0/pi,
     $        matmul(recip_lattice_conv,
     $        (kvec+dble((/ia, ib, ic/)))
     $        ))
            do ii=1, 3
              dblearray(ii,ia,ib,ic)=modulo(tempdble3(ii), 1.0d0)
            enddo
          enddo
        enddo
      enddo
      
      do ia=0, num_cell-1
        do ib=0, num_cell-1
          do ic=0, num_cell-1
            tempcnt=0
            do ja=0, num_cell-1
              do jb=0, num_cell-1
                do jc=0, num_cell-1
                  if (tempcnt.eq.0) then
                    do ii=1, 3
                      tempdble3(ii)=
     $                  modulo(dblearray(ii,ia,ib,ic)
     $                  -dblearray(ii,ja,jb,jc), 1.0d0)
                      if (tempdble3(ii) .gt. 0.5d0) then
                        tempdble3(ii)=tempdble3(ii)-1.0d0
                      endif
                    enddo
                    
c$$$                    write(179,'(7i5,9f12.6)')
c$$$     $                ia, ib,ic, ja, jb, jc,
c$$$     $                tempcnt, dblearray(:,ia,ib,ic),
c$$$     $                dblearray(:,ja,jb,jc),tempdble3
                    
                    if (dsqrt(sum(tempdble3**2)).lt. 1.0d-4)
     $                then
                      issame(:,ia,ib,ic)=(/ja,jb,jc/)
                      degennum(ja,jb,jc)=degennum(ja,jb,jc)+1
                      tempcnt=tempcnt+1
                    endif
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c$$$      write(179,*) 'dblearray'
c$$$      do ia=0, num_cell-1
c$$$        do ib=0, num_cell-1
c$$$          do ic=0, num_cell-1
c$$$            write(179,'(6f12.6,3i5)')
c$$$     $        kvec,dblearray(:,ia,ib,ic),issame(:,ia,ib,ic)
c$$$          enddo
c$$$        enddo
c$$$      enddo
      
      notzero=0
      do ia=0, num_cell-1
        do ib=0, num_cell-1
          do ic=0, num_cell-1
            if (degennum(ia,ib,ic) .ne. 0) then
              notzero=notzero+1
            endif
          enddo
        enddo
      enddo
c$$$      write(179,*) 'notzero', notzero
      if (
     $  (notzero .ne. num_cell)
     $  .or.
     $  (sum(degennum) .ne. num_cell**3)
     $  )
     $  then
c$$$        write(179,*) "zone folding wrong"
        CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
        CALL MPI_FINALIZE(mpierr)        
      else
        tempcnt=0
        do ia=0, num_cell-1
          do ib=0, num_cell-1
            do ic=0, num_cell-1            
              if (
     $          (issame(1,ia,ib,ic) .eq. ia) .and.
     $          (issame(2,ia,ib,ic) .eq. ib) .and.
     $          (issame(3,ia,ib,ic) .eq. ic)
     $          ) then
                tempcnt=tempcnt+1
                gvec(:,tempcnt)=(/ia,ib,ic/)
              endif
            enddo
          enddo
        enddo
      endif
      end


      subroutine hermitianeigen_cmplxdouble(dimen, w, datamat)
      implicit none
      
      integer, intent(in) :: dimen
      complex*16, intent(inout) :: datamat(dimen, dimen)
      double precision, intent(out) :: w(dimen)
      integer :: lwork, info, j, i, ifail, lwkopt
      complex*16, allocatable :: work(:)
      complex*16 :: worktemp(2*dimen)
      double precision :: rwork(dimen*3), rcondz(dimen), zerrbd(dimen),
     $  eerrbd, eps
      
      call zheev('V','U',dimen,datamat,dimen, w,worktemp,-1,rwork, info)
      lwkopt = worktemp(1)
      allocate(work(lwkopt))
      call zheev('V','U',dimen,datamat,dimen, w,work,lwkopt,rwork, info)
      deallocate(work)
      if (info .ne. 0) then
        write(*,*) 'error in diagonalization', info
        stop
      end if
      end subroutine hermitianeigen_cmplxdouble      


      subroutine check_if_unitary(ndim,ndim0, mat, fail_ind)
      implicit none
      integer, intent(in) :: ndim,ndim0
      complex*16, intent(in) :: mat(ndim0,ndim0)
      integer, intent(out) :: fail_ind

      integer :: ii,jj
      complex*16 :: mat2(ndim,ndim)
      

      mat2=0.0d0
      call zgemm('c','n',ndim,ndim,ndim,
     $  (1.0d0,0.0d0),mat,ndim0,mat,ndim0,
     $  (0.0d0,0.0d0),mat2,ndim)

      fail_ind=0
      do ii=1, ndim
        do jj=1, ndim
          if (ii .eq. jj) then
            if (cdabs(mat2(ii,jj)-1.0d0) .gt. 1.0d-8) then
              fail_ind=fail_ind+1
            endif
          else
            if (cdabs(mat2(ii,jj)) .gt. 1.0d-8) then
              fail_ind=fail_ind+1
            endif
          endif
        enddo
      enddo
      end


      double precision function find_fermi_level (a,b,machep,t,f,
     $  numk,ndim_k,n_mpi_k,ndim_k_max,kvec,nproc,nshift)

!*****************************************************************************80
!     
!     ! ZERO seeks the root of a function F(X) in an interval [A,B].
!     
!     Discussion:
!     
!     The interval [A,B] must be a change of sign interval for F.
!     That is, F(A) and F(B) must be of opposite signs.  Then
!     assuming that F is continuous implies the existence of at least
!     one value C between A and B for which F(C) = 0.
!     
!     The location of the zero is determined to within an accuracy
!     of 6 * MACHEPS * abs ( C ) + 2 * T.
!     
!     Thanks to Thomas Secretin for pointing out a transcription error in the
!     setting of the value of P, 11 February 2013.
!     
!     Licensing:
!     
!     This code is distributed under the GNU LGPL license. 
!     
!     Modified:
!     
!     11 February 2013
!     
!     Author:
!     
!     Original FORTRAN77 version by Richard Brent.
!     FORTRAN90 version by John Burkardt.
!     
!     Reference:
!     
!     Richard Brent,
!     Algorithms for Minimization Without Derivatives,
!     Dover, 2002,
!     ISBN: 0-486-41998-3,
!     LC: QA402.5.B74.
!     
!     Parameters:
!     
!     Input, real ( kind = 8 ) A, B, the endpoints of the change of 
!     sign interval.
!     
!     Input, real ( kind = 8 ) MACHEP, an estimate for the relative machine
!     precision.
!     
!     Input, real ( kind = 8 ) T, a positive error tolerance.
!     
!     Input, external real ( kind = 8 ) F, the name of a user-supplied
!     function, of the form "FUNCTION F ( X )", which evaluates the
!     function whose zero is being sought.
!     
!     Output, real ( kind = 8 ) ZERO, the estimated value of a zero of
!     the function F.
!     
      implicit none
      include "mpif.h"      

      integer, intent(in) :: numk, nproc,ndim_k(nproc), n_mpi_k(nproc),
     $  ndim_k_max
      double precision, intent(in) :: kvec(3,numk)
      
      double precision, intent(in) :: a,b,machep,t,nshift
      double precision, external :: f

      integer :: ii, maxcnt
      double precision ::  c,d,e,fa,fb,fc,m,p,q,r,s,sa,sb,tol
!     
!     Make local copies of A and B.
!     
      find_fermi_level=0.0d0
      sa = a
      sb = b
      fa = f (sa,numk,ndim_k,n_mpi_k,ndim_k_max,kvec,nshift)
      fb = f (sb,numk,ndim_k,n_mpi_k,ndim_k_max,kvec,nshift)

      c = sa
      fc = fa
      e = sb - sa
      d = e

      maxcnt=300
      do ii=1, maxcnt
        if ( abs ( fc ) < abs ( fb ) ) then
          sa = sb
          sb = c
          c = sa
          fa = fb
          fb = fc
          fc = fa
        end if
        tol = 2.0D+00 * machep * abs ( sb ) + t
        m = 0.5D+00 * ( c - sb )
        if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
          exit
        end if
        if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then
          e = m
          d = e
        else
          s = fb / fa
          if ( sa == c ) then
            p = 2.0D+00 * m * s
            q = 1.0D+00 - s
          else
            q = fa / fc
            r = fb / fc
            p = s * ( 2.0D+00 * m * q * ( q - r ) -
     $        ( sb - sa ) * ( r - 1.0D+00 ) )
            q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )
          end if
          if ( 0.0D+00 < p ) then
            q = - q
          else
            p = - p
          end if

          s = e
          e = d
          
          if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. 
     &      p < abs ( 0.5D+00 * s * q ) ) then
            d = p / q
          else
            e = m
            d = e
          end if
          
        end if
        
        sa = sb
        fa = fb
        
        if ( tol < abs ( d ) ) then
          sb = sb + d
        else if ( 0.0D+00 < m ) then
          sb = sb + tol
        else
          sb = sb - tol
        end if
        
        fb = f ( sb,numk,ndim_k,n_mpi_k,ndim_k_max,kvec,nshift)
        
        if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or.
     $    ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
          c = sa
          fc = fa
          e = sb - sa
          d = e
        end if
        
      end do
      
      if (ii .eq. maxcnt) then
        write(993,*) 'fermi level search fail'
        call ending                
      endif      
c$$$  write(*,*) 'maxcnt', ii
      
      find_fermi_level = sb
      
      end

! is_orthonormal_dcmplx calculates mata^\dagger matb and checks that it is equal to one.
      subroutine is_orthonormal_dcmplx(mata,lda,matb,ldb,matsize1,
     $  matsize2,matsize3,is_identity, string)

!     assumes that matsize3 < lda and ldb
!     multiplication of dconjg(transpose(a(1:matsize3,1:matsize1))) and b(1:matsize3,1:matsize2)
      
      implicit none
      
      integer,intent(in) :: matsize1,matsize2,matsize3,is_identity,
     $  lda,ldb
      complex*16,intent(in) :: mata(lda,matsize1),
     $  matb(ldb,matsize2)
      character(len=*), intent(in) :: string
      integer :: ii, jj
      double precision :: err, tempval
      complex*16 :: tempmat(matsize1,matsize2)

      call zgemm('c','n',matsize1,matsize2,matsize3,
     $  (1.d0,0.d0),mata,lda,
     &  matb,ldb,(0.d0,0.d0),
     $  tempmat,matsize1)

      if (is_identity .eq. 1) then
        tempval=1.0d0
      else
        tempval=0.0d0
      endif
      
      err=1.0d-6
      do ii=1, matsize1
        do jj=1, matsize2        
          if (ii .eq. jj) then
            if (cdabs(tempmat(ii,jj)-tempval) .gt. err) then
c$$$  if (maswrk) then
c$$$  write(993,'(a, 2i5,2f12.6)')
c$$$  $            'orthonormalization fail: '//trim(string),
c$$$  $            ii, jj, tempmat(ii,jj)
c$$$  call flush(993)
c$$$  endif
c$$$              write(179,'(a, 2i5,2f12.6)')
c$$$     $          'orthonormalization fail: '//trim(string),
c$$$     $          ii, jj, tempmat(ii,jj)
c$$$              call flush(179)              
              call ending
            endif
          else
            if (cdabs(tempmat(ii,jj)) .gt. err) then
c$$$  if (maswrk) then
c$$$  write(993,'(a, 2i5,2f12.6)')
c$$$  $            'orthonormalization fail: '//trim(string),
c$$$  $            ii, jj, tempmat(ii,jj)
c$$$  call flush(993)                              
c$$$  endif
c$$$              write(179,'(a, 2i5,2f12.6)')
c$$$     $          'orthonormalization fail: '//trim(string),
c$$$     $          ii, jj, tempmat(ii,jj)
c$$$              call flush(179)                            
              call ending              
            endif
          endif
        enddo
      enddo
      end



      subroutine check_hermitian(mata,matsize1,string)

!     assumes that matsize3 < lda and ldb
!     multiplication of dconjg(transpose(a(1:matsize3,1:matsize1))) and b(1:matsize3,1:matsize2)
      
      implicit none
      
      integer,intent(in) :: matsize1
      complex*16,intent(in) :: mata(matsize1,matsize1)
      character(len=*), intent(in) :: string
      integer :: ii, jj
      double precision :: err

      err=1.0d-6
      
      do ii=1, matsize1
        do jj=ii, matsize1        
          if (cdabs(mata(ii,jj)-dconjg(mata(jj,ii))) .gt. err) then
c$$$            write(179,'(a, 2i5,2f12.6, 2x, 2f12.6)')
c$$$     $        'hermitian fail: '//trim(string),
c$$$     $        ii, jj, mata(ii,jj),mata(jj,ii)
c$$$            call flush(179)
            write(993,'(a, 2i5,2f12.6, 2x, 2f12.6)')
     $        'hermitian fail: '//trim(string),
     $        ii, jj, mata(ii,jj),mata(jj,ii)
            call flush(993)            
            call ending
          endif
        enddo
      enddo
      end                                

! cal_mat_eig_co calculates h_nonint_co, and eigenvalues_co, and eigenvalues_co_all, and w_matrix_co. It would be trivial if zinv_m1_exist were false, i.e. if zinv_m1.dat did not exist, which is the case when doing dft+dmft.
! h_nonint_co is used only in cal_n_matrix, and w_matrix_co is used by write_den_matrix.
! eigenvalues_co=eigenvalues_co_all and w_matrix_co are the eigenvalues and eigenvectors of v_matrix^\dagger * eigenvalues * v_matrix
! h_nonint_co(k) = \sqrt{zinv}(k) * eigenvalues_co(k) * \sqrt{zinv}(k)
! zinv = 1 + \sum_i proj_renorm_co^\dagger(i,k) * zinv_m1(i) * proj_renorm_co(i,k)
! If zinv_m1.dat doesn't exist that means that zinv_m1 = 0 and h_nonint_co = eigenvalues_co.
! output: h_nonint_co
! output: eigenvalues_co, eigenvalues_co_all
! output: w_matrix_co
! input: v_matrix
! input: eigenvalues
! input: zinv_m1
! input: proj_renorm_co
      subroutine cal_mat_eig_co
      use comlowh_mod      
      implicit none
      include "mpif.h"      
      integer :: ibnd,iwan,jwan,ik,irk,iwancell,icell,iatom,ii,
     $  iorb,jorb,iio,jwancell,kwancell
      double precision :: std_dble,
     $  diag(num_wancell)
      complex*16 :: h_k_wan_co(num_wann,num_wann),tempvec(num_wann),
     $  mat_loc_ackw(num_orb_max,num_orb_max),
     $  tempmat1(num_orb_max,num_wancell),
     $  zinv(num_wancell,num_wancell),
     $  zinv_half(num_wancell,num_wancell)

      h_nonint_co=0.0d0
      do ik=1, ndim_kk_co(me1)
        irk=n_mpi_kk_co(me1)+ik

! h_k_wan_co = v_matrix^\dagger * eigenvalues * v_matrix
        h_k_wan_co=0.0d0
        do iwan=1, num_wann
          do jwan=1, num_wann
            do ibnd=1, num_bands                
              h_k_wan_co(iwan,jwan)
     $          =h_k_wan_co(iwan,jwan)
     $          +eigenvalues(ibnd,irk)
     $          *v_matrix(ibnd,jwan,irk)
     $          *dconjg(v_matrix(ibnd,iwan,irk))
            enddo
          enddo
        enddo

! eigenvalues_co_all=eigenvalues_co and h_k_wan_co are the eigenvalues and eigenvectors of v_matrix^\dagger * eigenvalues * v_matrix
        call hermitianeigen_cmplxdouble(num_wann,
     $    eigenvalues_co(1,ik),h_k_wan_co)

        eigenvalues_co_all(:,irk)=eigenvalues_co(:,ik)        
        w_matrix_co(:,:,ik)=h_k_wan_co

      enddo

      call mpi_allreduce_dble(eigenvalues_co_all,
     $  nqdiv*num_wancell,mpi_sum,mpi_comm_world)                    


      do ik=1, ndim_kk_co(me1)
        irk=n_mpi_kk_co(me1)+ik

        if (zinv_m1_exist) then
          zinv=0.0d0
          do iwancell=1, num_wancell
            zinv(iwancell,iwancell)=1.0d0
          enddo
          
          do icell=1, num_cell ! num_cell = 1
            do iatom=1, num_atom

              ii=iac(iatom,icell)
              mat_loc_ackw=0.0d0
              do iorb=1, num_orb(iatom)
                do jorb=1, num_orb(iatom)
                  iio=ioac(iorb,jorb,ii)
                  if (iio .ne. 0) then
                    mat_loc_ackw(iorb,jorb)
     $                =mat_loc_ackw(iorb,jorb)
     $                +zinv_m1(iio,ii)
                  endif
                enddo
              enddo
              
              tempmat1=0.0d0
              call zgemm('n','n',num_orb(iatom),num_wancell,
     $          num_orb(iatom),(1.d0,0.d0),
     $          mat_loc_ackw,num_orb_max,
     $          proj_renorm_co(1,1,iatom,icell,ik),num_orb_max,
     $          (0.d0,0.d0),
     $          tempmat1,num_orb_max)
              
              call zgemm('c','n',num_wancell,num_wancell,
     $          num_orb(iatom),(1.d0,0.d0),
     $          proj_renorm_co(1,1,iatom,icell,ik),num_orb_max,
     $          tempmat1,num_orb_max,(1.d0,0.d0),
     $          zinv,num_wancell)

            enddo ! iatom=1, num_atom
          enddo !  icell=1, num_cell ! num_cell = 1
          
          call hermitianeigen_cmplxdouble
     $      (num_wancell,diag,zinv)
          zinv_half=0.0d0
          do iwancell=1, num_wancell
            do jwancell=1, num_wancell
              do kwancell=1, num_wancell                              
                zinv_half(iwancell,jwancell)
     $            =zinv_half(iwancell,jwancell)
     $            +zinv(iwancell,kwancell)*dsqrt(diag(kwancell))
     $            *dconjg(zinv(jwancell,kwancell))
              enddo
            enddo
          enddo
          
          do iwancell=1, num_wancell
            do jwancell=1, num_wancell
              do kwancell=1, num_wancell                              
                h_nonint_co(iwancell,jwancell,ik)
     $            =h_nonint_co(iwancell,jwancell,ik)
     $            +zinv_half(iwancell,kwancell)
     $            *eigenvalues_co(kwancell,ik)
     $            *zinv_half(kwancell,jwancell)
              enddo
            enddo
          enddo
        else ! not (zinv_m1_exist)
          do iwancell=1, num_wancell
            h_nonint_co(iwancell,iwancell,ik)
     $        =eigenvalues_co(iwancell,ik)
          enddo
        endif
      enddo      

      end
      

! cal_hmat_rhomat_from_proj calculates and prints out a hamiltonian
!   and a density matrix. It does not return anything to its caller.
!  H = Nk^{-1} \sum_k P (E - mu) P^\dagger , rho = Nk^{-1} \sum_k P fermi(\beta(E - mu)) P^\dagger
! Nk is the number of k points, P = proj_in, E = eigval, mu is the Fermi level, fermi is the fermi distribution, beta is the inverse temperature.
      subroutine cal_hmat_rhomat_from_proj
     $  (numk,ndim_k,n_mpi_k,ndim_k_max,kvec,eigval,
     $  proj_in,textmessage)

      use comlowh_mod

      implicit none
      include "mpif.h"
      
      integer, intent(in) :: numk, ndim_k(nproc), n_mpi_k(nproc),
     $  ndim_k_max
      double precision, intent(in) :: kvec(3,numk),
     $  eigval(num_wancell,ndim_k_max)

      complex*16, intent(in) :: proj_in(num_orb_max,num_wancell,
     $  num_atom,num_cell,ndim_k_max)
      character(len=*), intent(in) :: textmessage            

      integer :: 
     $  iatom, isort,icell,
     $  atom_count,iorb,jorb,ik,irk,iband, ii,jj, k1,k2,k3,kk
      complex*16 :: 
     $  rhomat(num_orb_max,num_orb_max,num_atom,num_cell),
     $  hmat(num_orb_max,num_orb_max,num_atom,num_cell)      

      
      hmat=0.0d0
      rhomat=0.0d0
      
      do ik=1, ndim_k(me1)
        irk=n_mpi_k(me1)+ik
        do icell=1, num_cell
          do iatom=1, num_atom
            do iorb=1, num_orb(iatom)
              do jorb=1, num_orb(iatom)
                do iband=1, num_wancell
c$$$  
                  hmat(iorb,jorb,iatom,icell)
     $              =hmat(iorb,jorb,iatom,icell)+
     $              (eigval(iband,ik)-ef)
     $              *proj_in(iorb,iband,iatom,icell,ik)
     $              *dconjg(proj_in(jorb,iband,iatom,icell,ik))
     $              /dble(numk)

                  rhomat(iorb,jorb,iatom,icell)
     $              =rhomat(iorb,jorb,iatom,icell)+
     $              proj_in(iorb,iband,iatom,icell,ik)
     $              /(dexp(beta*(eigval(iband,ik)-ef))+1.0d0)
     $              *dconjg(proj_in(jorb,iband,iatom,icell,ik))
     $              /dble(numk)                  
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      
      call mpi_allreduce_dcmplx(rhomat,num_orb_max**2*num_atom*num_cell,
     $  mpi_sum,mpi_comm_world)
      call mpi_allreduce_dcmplx(hmat,num_orb_max**2*num_atom*num_cell,
     $  mpi_sum,mpi_comm_world)
      
      if (maswrk) then
        write(iun,'(a, i5, a, i5)')
     $    'hmat and rhomat from '//textmessage
        do icell=1, num_cell
          do iatom=1,num_atom
            write(iun,'(a, i5, a, i5)') 'rhomat, iatom  ',
     $        iatom, '  icell  ', icell
            do ii=1, num_orb(iatom)
              write(iun,'(100(2f12.6, 5x))')
     $          (rhomat(ii,iorb,iatom,icell), iorb=1,num_orb(iatom))
            enddo
          enddo
        enddo
        do icell=1, num_cell
          do iatom=1,num_atom
            write(iun,'(a, i5, a, i5)') 'hmat, iatom  ',
     $        iatom, '  icell  ', icell
            do ii=1, num_orb(iatom)
              write(iun,'(100(2f12.6, 5x))')
     $          (hmat(ii,iorb,iatom,icell), iorb=1,num_orb(iatom))
            enddo
          enddo
        enddo
      endif      
      end





      subroutine find_proc_ind(ndim_proc, n_mpi_proc, nproc,indin,
     $  proc_num,indout)
!     indin starting from 1
!     proc_num starting from 0
!     indout starting from 1      

      implicit none
      integer, intent(in) :: nproc, ndim_proc(nproc),n_mpi_proc(nproc),
     $  indin
      integer, intent(out) :: proc_num,indout
      integer :: ip
      
      proc_num=0
      indout=0
      do ip=0, nproc-1
        if (indin .ge. n_mpi_proc(ip+1)+1) then
          proc_num=ip
        endif
      enddo

      indout=indin-n_mpi_proc(proc_num+1)
      end

      subroutine write_basicarray(maxdim2,nindo)
      use comlowh_mod
!     cidim and nindo is different only if nl is different from 1
      implicit none
      integer, intent(out) :: maxdim2,nindo(num_atom)
      integer :: maxdim,maxsize,lval,iatom,ii,jj,ip,iq,
     $  cidim(num_atom),nl(num_atom),iso,
     $  ll_ind(num_atom),iorbital(num_atom),nind(num_atom),
     $  csize(num_atom),cix(num_atom),ip1

      double precision :: vol,rbas(3,3),dbletemp(3)
      integer, allocatable :: isx(:,:),sigind(:,:,:)
      

      rbas=real_lattice/0.5291772106712d0
      call dcross(dbletemp,rbas(:,2),rbas(:,3))
      vol=sum(dbletemp*rbas(:,1))
      nindo=0
      if (.not. maswrk) return
      
      call set_basic_array_dimension(maxdim,maxdim2)
      
      allocate(isx(maxdim,num_atom))
      isx=0
      allocate(sigind(maxdim,maxdim,num_atom))
      sigind=0
      
      maxsize=maxval(ioac)      
      
      do iatom=1, num_atom
        lval=alimj(2,1,wan_correlated(1,iatom))
        if (is_spinorbit .eq. 0) then
          nindo(iatom)=2*lval+1
          cidim(iatom)=2*lval+1
          nind(iatom)=2*lval+1          
        else
          nindo(iatom)=(2*lval+1)*2
          cidim(iatom)=(2*lval+1)*2
          nind(iatom)=(2*lval+1)*2
        endif
        
        nl(iatom)=1
        ll_ind(iatom)=lval
        iorbital(iatom)=iatom
        cix(iatom)=iatom
        csize(iatom)=maxval(ioac(:,:,iac(iatom,1)))
        do ii=1, maxdim
          isx(ii,iatom)=ii
        enddo

        do ii=1, maxdim
          do jj=1, maxdim
            sigind(:,:,iatom)=ioac(:,:,iac(iatom,1))
          enddo
        enddo
      enddo
      
      open(unit=990, file='BasicArrays.dat')
      WRITE(990, *) 'nat, iso, norbitals, ncix, natom'
      if (is_spinorbit .eq. 1) then
        iso=2
      else
        iso=1
      endif
      WRITE(990, *) 1,iso,num_atom,num_atom,num_atom
      WRITE(990, *) 'nkpt, nmat, nume'
      WRITE(990, *) nqdiv_fi, num_wann, num_wann
      WRITE(990, *) 'Qcomplex'
      WRITE(990, *) 'True'
      WRITE(990, *) 'lmax2, maxdim2, maxdim, maxsize'
      WRITE(990, *) 4, maxdim2,maxdim,maxsize
      WRITE(990, *) 'nindo'
      WRITE(990, *) (nindo(iatom), iatom=1, num_atom)
      WRITE(990, *) 'cidim'
      WRITE(990, *) (cidim(iatom), iatom=1, num_atom)
      WRITE(990, *) 'nl'
      WRITE(990, *) (nl(iatom), iatom=1, num_atom)
      WRITE(990, *) 'll, iorbital, cix, nind'
      do iatom=1,num_atom
        WRITE(990, *) ll_ind(iatom)
        WRITE(990, *) iorbital(iatom)
        WRITE(990, *) cix(iatom)
        WRITE(990, *) nind(iatom)
      enddo
      WRITE(990, *) 'csize'
      WRITE(990, *) (csize(iatom), iatom=1,num_atom) ! independent number
      
      WRITE(990, *) 'iSx'
      do iatom=1,num_atom
        do ip1=1,nindo(iatom)
          WRITE(990, *) iSx(ip1,iatom) ! nonzero block number
        enddo
      enddo
      WRITE(990, *) 'Sigind'
      DO iatom=1,num_atom
        DO ip=1,cidim(iatom)
          DO iq=1,cidim(iatom)
            WRITE(990, *) Sigind(ip,iq,iatom)
          ENDDO
        ENDDO
      ENDDO
      
      WRITE(990, *) 'EF, VOL'
      WRITE(990, *) ef/13.60569300984d0, VOL ! rydberg,volume in bohr
      close(990)
      
      deallocate(isx)
      deallocate(sigind)
      end


      subroutine set_basic_array_dimension(maxdim,maxdim2)
      
      use comlowh_mod

      implicit none
      integer, intent(out) :: maxdim,maxdim2

      integer :: iatom,lval
      
      maxdim=0
      maxdim2=0
      
      do iatom=1, num_atom
        lval=alimj(2,1,wan_correlated(1,iatom))
        if (is_spinorbit .eq. 0) then
          maxdim2=maxval((/maxdim2,2*lval+1/))
        else
          maxdim2=maxval((/maxdim2,4*lval+2/))
        endif
      enddo
      maxdim=maxdim2
      end
      


      subroutine prepare_optics
      
      use comlowh_mod
      implicit none
      include "mpif.h"      

      integer :: ik,irk,icell,iatom,ind,ibnd,jbnd,proc_num,indout,
     $  nindo(num_atom),maxdim2,ii,iio,iomega,idir,iwan,jwan,
     $  i1,j1,tmpint(6)
      integer :: irequest,istatus(mpi_status_size)     
      logical :: mpiflag
      double precision :: eigenvalues_fi_all(num_wancell,nqdiv_fi),
     $  tempdble1,tempdble2,pi
      complex*16 ::ai,
     $  tempmat(num_wancell),
     $  tempmat2(num_orb_max,num_wancell,num_atom,num_cell),
     $  momentum_wan(num_wancell,num_wancell,nqdiv,3),
     $  momentum_fine(num_wancell,num_wancell,ndim_kk_fi_max,3),
     $  mat_r_wan(nrpts),
     $  tempmat1(num_wann,num_wann),tempmat3(num_wann,num_wann,3)
      character*30 :: ch_atom
      
!     basicarrays

      pi=4.0d0*datan2(1.0d0,1.0d0)
      ai=dcmplx(0.0d0, 1.0d0)      
      call write_basicarray(maxdim2,nindo)      
      
!     sig.inp*
      do iatom=1, num_atom
        ii=iac(iatom,1)
        write(ch_atom,*) iatom
        open(unit=56, file='sig.inp'//trim(adjustl(ch_atom)))

        do iomega=1, n_omega
          write(56,'(f20.12,4x)', advance='no') frequency_dble(iomega)
          do iio=1, maxval(ioac(:,:,ii))
            write(56,'(2f20.12,4x)',advance='no')
     $        sig(iio,ii,iomega)-dc(iio,ii)
          enddo
          write(56,*)
        enddo
        close(56)
      enddo
      

!     udmft0
      
      tempmat=0.0d0
      
      
      if (.not. maswrk) then
        do ik=1, ndim_kk_fi(me1)
          irk=n_mpi_kk_fi(me1)+ik
          call mpi_isend(proj_renorm(1,1,1,1,ik),
     $      num_orb_max*num_wancell*num_atom*num_cell,
     $      mpi_double_complex,0,irk,
     $      mpi_comm_world,irequest,mpierr)
          call mpi_wait(irequest,istatus,mpierr)
          call mpi_test(irequest,mpiflag,istatus,mpierr)
          if (.not.mpiflag)
     $      call mpi_abort(mpi_comm_world,mpierr)
        enddo
      else
        open(unit=990, file='Udmft.0',form='unformatted')
        open(unit=991, file='Udmft.0_f')        
        write(990) nqdiv_fi,1,num_atom
        write(991,'(3i5)') nqdiv_fi,1,num_atom        
        write(990) (nindo(iatom), iatom=1,num_atom)
        write(991,'(1000i5)') (nindo(iatom), iatom=1,num_atom)        
        do irk=1, nqdiv_fi
          call find_proc_ind(ndim_kk_fi, n_mpi_kk_fi, nproc,irk,
     $      proc_num,indout)
          if (proc_num .eq. 0)then
            
            write(990) irk, num_wann, maxdim2, num_atom, 1
            write(991,'(5i5)') irk, num_wann, maxdim2, num_atom, 1            
            write(990) 1
            write(991,'(i5)') 1            
            do icell=1, num_cell
              DO iatom=1,num_atom ! atom
                DO ind=1,nindo(iatom)
                  if (ind .le. num_orb(iatom)) then
                    write(990)
     $                (dconjg(proj_renorm(ind,ibnd,iatom,icell,indout)),
     $                ibnd=1,num_wancell)
                    write(991, '(10000(2f12.6,2x))')
     $                (dconjg(proj_renorm(ind,ibnd,iatom,icell,indout)),
     $                ibnd=1,num_wancell)                    
                  else
                    write(990)
     $                (tempmat(ibnd),ibnd=1,num_wancell)
                    write(991, '(10000(2f12.6,2x))')
     $                (tempmat(ibnd),ibnd=1,num_wancell)                    
                  endif
                enddo
              ENDDO
            ENDDO
          else
            tempmat2=0.0d0
            call mpi_irecv(tempmat2,
     $        num_orb_max*num_wancell*num_atom*num_cell,            
     $        mpi_double_complex,proc_num,
     $        irk,mpi_comm_world,irequest,mpierr)
            call mpi_wait(irequest,istatus,mpierr)
            call mpi_test(irequest,mpiflag,istatus,mpierr)
            if (.not.mpiflag)
     $        call mpi_abort(mpi_comm_world,mpierr)
            write(990) irk, num_wann, maxdim2, num_atom, 1
            write(991,'(5i5)') irk, num_wann, maxdim2, num_atom, 1                        
            write(990) 1
            write(991,'(i5)') 1                        
            do icell=1, num_cell
              DO iatom=1,num_atom ! atom
                DO ind=1,nindo(iatom)
                  if (ind .le. num_orb(iatom)) then
                    write(990)
     $                (dconjg(tempmat2(ind,ibnd,iatom,icell)),
     $                ibnd=1,num_wancell)
                    write(991, '(10000(2f12.6,2x))')
     $                (dconjg(tempmat2(ind,ibnd,iatom,icell)),
     $                ibnd=1,num_wancell)                    
                  else
                    write(990)
     $                (tempmat(ibnd),ibnd=1,num_wancell)
                    write(991, '(10000(2f12.6,2x))')
     $                (tempmat(ibnd),ibnd=1,num_wancell)                    
                  endif
                  
                enddo
              ENDDO
            ENDDO
          endif
        enddo
        close(990)
        close(991)
      endif

      call mpi_barrier(mpi_comm_world,mpierr)

!     eigenvalues
      eigenvalues_fi_all=0.0d0
      do ik=1, ndim_kk_fi(me1)
        irk=n_mpi_kk_fi(me1)+ik        
        eigenvalues_fi_all(:,irk)=eigenvalues_fi(:,ik)
      enddo
      
      call mpi_allreduce_dble
     $  (eigenvalues_fi_all,nqdiv_fi*num_wann,mpi_sum,mpi_comm_world)               
      if (maswrk) then
c$$$  case.energy        
        open(unit=113, file=adjustl(trim(seed_name))//'_optics.energy')
        write(113, *)
        write(113, *)
        do ik=1, nqdiv_fi
          write(113,'(3e19.12,i10,2i6,e19.12)')
     $      (kpt_latt_fi(ii,ik), ii=1, 3),
     $      ik, num_wann,num_wann, 1.0d0
          do ibnd=1, num_wann
            write(113, *) ibnd, eigenvalues_fi_all(ibnd,ik)
     $        /13.60569300984d0
          enddo
        enddo
        close(113)

c$$$  case.symop
        open(unit=113, file=adjustl(trim(seed_name))//'_optics.symop')
        write(113, '(i6)')   1
        write(113, '(3f8.5)')   1.0d0, 0.0d0, 0.0d0
        write(113, '(3f8.5)')   0.0d0, 1.0d0, 0.0d0
        write(113, '(3f8.5)')   0.0d0, 0.0d0, 1.0d0
        write(113,*)
        write(113,*)         
        close(113)
      endif
c$$$  case.mommat

      momentum_wan=0.0d0
      if (maswrk) then
        open(unit=55, file='momentum_optics_wan.dat')
        do idir=1,3
          do irk=1, nqdiv            
            do iwan=1, num_wann
              do jwan=1, num_wann          
                read(55,*) (tmpint(ii),ii=1,6),
     $            tempdble1,tempdble2
                momentum_wan(iwan,jwan,irk,idir)
     $            =dcmplx(tempdble1,tempdble2) 
              enddo
            enddo
          enddo
        enddo
        close(55)
      endif

      call mpi_bcast(momentum_wan,num_wann**2*nqdiv*3,
     $  mpi_double_complex,0,mpi_comm_world,mpierr)

      
      momentum_fine=0.0d0      
      do idir=1,3
        
        do iwan=1, num_wann
          do jwan=1, num_wann
            
            mat_r_wan=0.0d0
            
            do i1=1, nrpts
              do j1=1, nqdiv
                mat_r_wan(i1)=mat_r_wan(i1)
     $            +momentum_wan(iwan,jwan,j1,idir)*
     $            cdexp(-2.0*ai*pi*
     $            sum(irvec(:,i1)*kpt_latt(:,j1)))
     $            /dble(nqdiv)
              enddo
            enddo
c$$$  write(179,*) iwan,jwan, 'step2'
c$$$  call flush(179)          
            
            do ik=1, ndim_kk_fi(me1)
              irk=n_mpi_kk_fi(me1)+ik
              do j1=1, nrpts
                momentum_fine(iwan,jwan,ik,idir)
     $            =momentum_fine(iwan,jwan,ik,idir)
     $            +mat_r_wan(j1)*
     $            cdexp(2.0*ai*pi*
     $            sum(kpt_latt_fi(:,irk)
     $            *irvec(:,j1))
     $            )/dble(ndegen(j1))
              enddo
            enddo
c$$$  write(179,*) iwan,jwan, 'step3'
c$$$  call flush(179)          
          enddo
        enddo
      enddo
      

      if (.not. maswrk) then        
        do ik=1, ndim_kk_fi(me+1)
          irk=n_mpi_kk_fi(me+1)+ik
          tempmat3=0.0d0          
          do idir=1, 3            
            tempmat1=0.0d0
            call zgemm('n','n',num_wann,num_wann,num_wann,
     $        (1.d0,0.d0),momentum_fine(1,1,ik,idir),num_wann,
     $        w_matrix_fi(1,1,ik),num_wann,
     &        (0.d0,0.d0),tempmat1,num_wann)
            call zgemm('c','n',num_wann,num_wann,num_wann,
     $        (1.d0,0.d0),w_matrix_fi(1,1,ik),num_wann,
     $        tempmat1,num_wann,
     &        (0.d0,0.d0),tempmat3(1,1,idir),num_wann)
          enddo

          call mpi_isend(tempmat3,
     $      num_wann**2*3,mpi_double_complex,0,irk,
     $      mpi_comm_world,irequest,mpierr)
          call mpi_wait(irequest,istatus,mpierr)
          call mpi_test(irequest,mpiflag,istatus,mpierr)
          if (.not.mpiflag)
     $      call mpi_abort(mpi_comm_world,mpierr)
        enddo
      else
        open(unit=990, file=trim(adjustl(seed_name))
     $    //'_optics.mommat')
        open(unit=991, file=trim(adjustl(seed_name))
     $    //'_optics.mommat1')
        write(990, *) "  Complex momentum matrix elements: "//
     $    "(Re_x,Im_x),(Re_y,Im_y),(Re_z,Im_z)"
        write(991, *) "  Complex momentum matrix elements: "//
     $    "(Re_x,Im_x),(Re_y,Im_y),(Re_z,Im_z)"                  
        do irk=1, nqdiv_fi
          call find_proc_ind(ndim_kk_fi,n_mpi_kk_fi,nproc,
     $      irk,proc_num,indout)
          if (proc_num .eq. 0)then
            tempmat3=0.0d0          
            do idir=1, 3            
              tempmat1=0.0d0
              call zgemm('n','n',num_wann,num_wann,num_wann,
     $          (1.d0,0.d0),momentum_fine(1,1,indout,idir),num_wann,
     $          w_matrix_fi(1,1,indout),num_wann,
     &          (0.d0,0.d0),tempmat1,num_wann)
              call zgemm('c','n',num_wann,num_wann,num_wann,
     $          (1.d0,0.d0),w_matrix_fi(1,1,indout),num_wann,
     $          tempmat1,num_wann,
     &          (0.d0,0.d0),tempmat3(1,1,idir),num_wann)
            enddo
          else
            tempmat3=0.0d0
            call mpi_irecv(tempmat3,
     $        num_wann**2*3,mpi_double_complex,proc_num,
     $        irk,mpi_comm_world,irequest,mpierr)
            call mpi_wait(irequest,istatus,mpierr)
            call mpi_test(irequest,mpiflag,istatus,mpierr)
            if (.not.mpiflag)
     $        call mpi_abort(mpi_comm_world,mpierr)
          endif
          write(990, *)
          write(990, '(2x, a4,i6,a15, 2i5,a4,2f5.2,a3,i10)')
     $      ' KP:',irk,' NEMIN NEMAX : ',1, num_wann,' dE:',
     $      0.0d0, 0.0d0,' K:',  irk
          write(990, *)

          write(991, *)
          write(991, '(2x, a4,i6,a15, 2i5,a4,2f5.2,a3,i10)')
     $      ' KP:',irk,' NEMIN NEMAX : ',1, num_wann,' dE:',
     $      0.0d0, 0.0d0,' K:',  irk
          write(991, *)
          do ibnd=1, num_wann
            do jbnd=ibnd, num_wann
              write(990,'(3X,2I4,6E13.6,F13.8)')
     $          ibnd,jbnd,(tempmat3(ibnd,jbnd,idir), idir=1, 3),
     $          0.0d0

              write(991,'(3X,2I4,6E13.6,F13.8)')
     $          ibnd,jbnd,1.0d0,0.0d0,1.0d0,0.0d0,1.0d0,0.0d0,
     $          0.0d0              
            enddo
          enddo              
        enddo
        close(990)
        close(991)        
      endif

      call mpi_barrier(mpi_comm_world,mpierr)      
      end

! fermi_level_spectral_function_local calculates and prints out local_spectral_matrix_ef.dat, which is under certain circumstances used by comdmft.py. It also writes out n_loc.dat, which are never used by ComDMFT.
      subroutine fermi_level_spectral_function_local
      use comlowh_mod
      
      implicit none
      integer :: ii,iorb,jorb,iomega

      double precision :: pi
      complex*16 :: one_matrix_f(num_orb_max,num_orb_max),
     $  g2mat(num_orb_max, num_orb_max, n_iac_mat_i:n_iac_mat_f),
     $  g3mat(num_orb_max, num_orb_max, n_iac_mat_i:n_iac_mat_f),
     $  spectral(num_orb_max, num_orb_max, n_iac_mat_i:n_iac_mat_f),
     $  n_loc(num_orb_max, num_orb_max, n_iac_mat_i:n_iac_mat_f),
     $  deltaij

      pi=4.0d0*datan2(1.0d0,1.0d0)

      g2mat=0.0d0
      g3mat=0.0d0
      spectral=0.0d0
      n_loc=0.0d0
      
      do ii=n_iac_mat_i, n_iac_mat_f
        
        one_matrix_f=0.0d0
        do iorb=1, norb_iac(ii)
          one_matrix_f(iorb,iorb)=1.0d0
        enddo
        
        if ((is_magnetic .eq. 1) .and. (ii .eq. 0)) cycle

        g2mat(:,:,ii)=(gloc(:,:,n_omega,ii)
     $    +transpose(dconjg(gloc(:,:,n_omega,ii))))
     $    /2.0d0*frequency(n_omega)**2
        g3mat(:,:,ii)=(gloc(:,:,n_omega,ii)
     $    -transpose(dconjg(gloc(:,:,n_omega,ii)))
     $    -one_matrix_f*2.0/frequency(n_omega))
     $    /2.0d0*frequency(n_omega)**3                

        do iorb=1, norb_iac(ii)
          do jorb=1, norb_iac(ii)
            if (iorb .eq. jorb) then
              deltaij=1.0d0
            else
              deltaij=0.0d0
            endif
            do iomega=1, n_omega
              spectral(iorb,jorb,ii)
     $          =spectral(iorb,jorb,ii)
     $          +1.0d0/beta*cdexp(-beta/2.0d0*frequency(iomega))
     $          *(
     $          gloc(iorb,jorb,iomega,ii)
     $          -deltaij/frequency(iomega)
     $          -g2mat(iorb,jorb,ii)/(frequency(iomega))**2
     $          -g3mat(iorb,jorb,ii)/(frequency(iomega))**3                  
     $          )
     $          +1.0d0/beta*cdexp(beta/2.0d0*frequency(iomega))
     $          *(
     $          dconjg(gloc(jorb,iorb,iomega,ii))
     $          +deltaij/frequency(iomega)
     $          -g2mat(iorb,jorb,ii)/(frequency(iomega))**2
     $          +g3mat(iorb,jorb,ii)/(frequency(iomega))**3                 
     $          )

              n_loc(iorb,jorb,ii)
     $          =n_loc(iorb,jorb,ii)
     $          -1.0d0/beta*cdexp(-beta*frequency(iomega))
     $          *(
     $          gloc(iorb,jorb,iomega,ii)
     $          -deltaij/frequency(iomega)
     $          -g2mat(iorb,jorb,ii)/(frequency(iomega))**2
     $          -g3mat(iorb,jorb,ii)/(frequency(iomega))**3                  
     $          )
     $          -1.0d0/beta*cdexp(-beta*frequency(iomega))
     $          *(
     $          dconjg(gloc(jorb,iorb,iomega,ii))
     $          +deltaij/frequency(iomega)
     $          -g2mat(iorb,jorb,ii)/(frequency(iomega))**2
     $          +g3mat(iorb,jorb,ii)/(frequency(iomega))**3                 
     $          )              
            enddo
            spectral(iorb,jorb,ii)
     $        =spectral(iorb,jorb,ii)
     $        -deltaij/2.0d0
     $        -g3mat(iorb,jorb,ii)*beta**2/4.0d0
     $        *((1.0d0/2.0d0)**2-(1.0d0/2.0d0))

            n_loc(iorb,jorb,ii)
     $        =n_loc(iorb,jorb,ii)
     $        +deltaij/2.0d0
     $        -g2mat(iorb,jorb,ii)*beta/4.0d0
          enddo
        enddo
      enddo
      if (maswrk) then
! local_spectral_matrix_ef.dat is never used by ComDMFT.
        open(305, file='local_spectral_matrix_ef.dat')
        do ii=n_iac_mat_i, n_iac_mat_f
          if ((is_magnetic .eq. 1) .and. (ii .eq. 0)) cycle
          do iorb=1, norb_iac(ii)
            write(305, '(100(f20.12, 2x))')
     $        (dble(spectral(iorb,jorb,ii))
     $        /pi*beta*(-1.0d0), jorb=1, norb_iac(ii))
          enddo
        enddo
        close(305)
! n_loc.dat is never used by ComDMFT.
        open(305, file='n_loc.dat')
        do ii=n_iac_mat_i, n_iac_mat_f
          if ((is_magnetic .eq. 1) .and. (ii .eq. 0)) cycle
          do iorb=1, norb_iac(ii)
            write(305, '(100(2f20.12, 2x))')
     $        (n_loc(iorb,jorb,ii)
     $        , jorb=1, norb_iac(ii))
          enddo
        enddo
        close(305)
        
      endif
      end

! fermi_level_spectral_function computes and prints out  spectral_matrix_ef.dat
      subroutine fermi_level_spectral_function
     $  (numk, ndim_k, n_mpi_k,ndim_k_max, g_lat_kw_diagonal)
      
      use comlowh_mod
      
      implicit none
      include "mpif.h"      

      integer, intent(in) :: numk, ndim_k(nproc), n_mpi_k(nproc),
     $  ndim_k_max      

      complex*16, intent(in) ::
     $  g_lat_kw_diagonal(num_wancell, ndim_k_max, n_omega)
      integer :: iomega,ik,irk,ibnd

      double precision :: pi
      complex*16 :: g2, g3,spectral(num_wancell,numk)

      pi=4.0d0*datan2(1.0d0,1.0d0)

      g2=0.0d0
      g3=0.0d0
      spectral=0.0d0

      do ik=1, ndim_k(me1)
        irk=n_mpi_k(me1)+ik
        do ibnd=1, num_wancell
          
          g2=dble(g_lat_kw_diagonal(ibnd,ik,n_omega))
     $      *frequency(n_omega)**2
          g3=(dimag(g_lat_kw_diagonal(ibnd,ik,n_omega))
     $      -1.0d0/frequency(n_omega))
     $      *frequency(n_omega)**3                
          
          do iomega=1, n_omega
            spectral(ibnd,irk)
     $        =spectral(ibnd,irk)
     $        +1.0d0/beta*cdexp(-beta/2.0d0*frequency(iomega))
     $        *(
     $        g_lat_kw_diagonal(ibnd,ik,iomega)
     $        -1.0d0/frequency(iomega)
     $        -g2/(frequency(iomega))**2
     $        -g3/(frequency(iomega))**3                  
     $        )
     $        +1.0d0/beta*cdexp(beta/2.0d0*frequency(iomega))
     $        *(
     $        dconjg(g_lat_kw_diagonal(ibnd,ik,iomega))
     $        +1.0d0/frequency(iomega)
     $        -g2/(frequency(iomega))**2
     $        +g3/(frequency(iomega))**3                 
     $        )
          enddo
          spectral(ibnd,irk)
     $      =spectral(ibnd,irk)
     $      -1.0d0/2.0d0
     $      -g3*beta**2/4.0d0
     $      *((1.0d0/2.0d0)**2-(1.0d0/2.0d0))
        enddo
      enddo

      call mpi_allreduce_dcmplx(spectral,
     $  num_wancell*numk,mpi_sum,mpi_comm_world)      
      
      if (maswrk) then
        open(305, file='spectral_matrix_ef.dat')
        do ibnd=1, num_wancell
          write(305, '(100000(f20.12, 2x))')
     $      (dble(spectral(ibnd,irk))
     $      /pi*beta*(-1.0d0), irk=1, numk)
        enddo
        close(305)
      endif
      end


! fermi_level_spectral_function_ksum computes and prints out  spectral_matrix_ef.dat
      subroutine fermi_level_spectral_function_ksum
     $  (numk, ndim_k, n_mpi_k,ndim_k_max, g_lat_kw_diagonal)
      
      use comlowh_mod
      
      implicit none
      include "mpif.h"      

      integer, intent(in) :: numk, ndim_k(nproc), n_mpi_k(nproc),
     $  ndim_k_max      

      complex*16, intent(in) ::
     $  g_lat_kw_diagonal(num_wancell, ndim_k_max, n_omega)
      integer :: iomega,ik,irk,ibnd

      double precision :: pi
      complex*16 :: g2, g3,spectral(num_wancell),
     $  g_lat_w(num_wancell,n_omega)

      pi=4.0d0*datan2(1.0d0,1.0d0)

      g2=0.0d0
      g3=0.0d0
      spectral=0.0d0


      g_lat_w=0.0d0
      do ik=1, ndim_k(me1)
        do iomega=1, n_omega
          do ibnd=1, num_wancell
            g_lat_w(ibnd,iomega)=g_lat_w(ibnd,iomega)
     $        +g_lat_kw_diagonal(ibnd,ik,iomega)/numk
          enddo
        enddo
      enddo

      call mpi_allreduce_dcmplx(g_lat_w,
     $  num_wancell*n_omega,mpi_sum,mpi_comm_world)      
      
            
      do ibnd=1, num_wancell
        
        g2=dble(g_lat_w(ibnd,n_omega))
     $    *frequency(n_omega)**2
        g3=(dimag(g_lat_w(ibnd,n_omega))
     $    -1.0d0/frequency(n_omega))
     $    *frequency(n_omega)**3                
        
        do iomega=1, n_omega
          spectral(ibnd)
     $      =spectral(ibnd)
     $      +1.0d0/beta*cdexp(-beta/2.0d0*frequency(iomega))
     $      *(
     $      g_lat_w(ibnd,iomega)
     $      -1.0d0/frequency(iomega)
     $      -g2/(frequency(iomega))**2
     $      -g3/(frequency(iomega))**3                  
     $      )
     $      +1.0d0/beta*cdexp(beta/2.0d0*frequency(iomega))
     $      *(
     $      dconjg(g_lat_w(ibnd,iomega))
     $      +1.0d0/frequency(iomega)
     $      -g2/(frequency(iomega))**2
     $      +g3/(frequency(iomega))**3                 
     $      )
        enddo
        spectral(ibnd)
     $    =spectral(ibnd)
     $    -1.0d0/2.0d0
     $    -g3*beta**2/4.0d0
     $    *((1.0d0/2.0d0)**2-(1.0d0/2.0d0))
      enddo

      
      if (maswrk) then
        open(305, file='spectral_matrix_ef.dat')
        do ibnd=1, num_wancell
          write(305, '(f20.12, 2x)')
     $      dble(spectral(ibnd))
     $      /pi*beta*(-1.0d0)
        enddo
        close(305)
      endif
      end      

      subroutine indexing(ntot,ndivision,divisionarray,flag,n1,n2)
      use comlowh_mod
      implicit none
      integer, intent(in) :: ntot,ndivision,
     $  divisionarray(ndivision),flag
      integer, intent(inout) :: n1, n2(ndivision)
      integer :: ii,jj,tempcnt,tmpsize

      tmpsize=1
      do ii=1, ndivision
        tmpsize=tmpsize*divisionarray(ii)
      enddo
      if (tmpsize .ne.  ntot) then
        if (maswrk) then
          write(993,*) 'array_division wrong'
          write(993,'(100(i10,2x))') tmpsize, ntot,
     $    divisionarray
          call flush(993)
        endif
c$$$        write(179,*) 'array_division wrong'        
c$$$        write(179,'(100(i10,2x))') tmpsize, ntot,
c$$$     $    divisionarray
        call ending
      endif
      if (flag .eq. 1) then
        n1=n2(1)
        do ii=2, ndivision
          tempcnt=1
          do jj=1, ii-1
            tempcnt=tempcnt*divisionarray(jj)
          enddo
          n1=n1+(n2(ii)-1)*tempcnt
        enddo
      else
        n2=0
        tempcnt=n1
        do ii=1, ndivision-1
          n2(ii)=tempcnt-((tempcnt-1)/divisionarray(ii))
     $      *divisionarray(ii)
          tempcnt=(tempcnt-n2(ii))/divisionarray(ii)+1
        enddo
        n2(ndivision)=tempcnt
      endif
      end

! cal_n_matrix is almost exactly the same as cal_n_lattice+cal_tos_analytical_tail_fun, except 
!   that the former is evaluated on the coarse k-point grid associated with ndiv and nqdiv,
!   while the latter is evaluated on the fine k-point grid associated with ndiv_fi and nqdiv_fi.
!   The former is used to produce the density matrix which is supplied to flapwmbpt, while
!   the latter is used to produce an electron count for adjusting the Fermi level.
! cal_n_matrix calculates the n_matrix:
! It calculates g_lat_inv_co(omega) = frequency(iomega) - h_nonint_co(:,:,ik) + \sum_{atom} proj_renorm_co(k)^\dagger * (sig(omega) - dc - I * gamma_selfenergy) * proj_renorm_co(k)
! Next it calculates gmat = (ef + g_lat_inv_co(:,:,n_omega))^{-1}
! Next it uses gmat   to calculate g2mat and g3mat, which will help with the tail at large momenta.
! Next, it loops over omega, and in this loop it recalculates gmat = (ef + g_lat_inv_co(:,:,omega))^{-1}, and then uses gmat, g2mat, g3mat to add a contribution to n_matrix.
! All of this logic has updated only the P_wann sector of n_matrix. Next we update the P_rest sector of of n_matrix, by adding   fermi(beta(eigenvalues_co_rest(k) -ef)), where eigenvalues_co_rest is a diagonal matrix and is defined only in the P_rest sector.
! Last it repeats the above in a loop over k, so n_matrix is summed over k.
! ---------------
! output: n_matrix(num_bands,num_bands,nqdiv)
! local: g_lat_inv_co(num_wann,num_wann,n_omega)
! input: h_nonint_co(num_wann,num_wann,ndim_kk_co_max)
! input: proj_renorm_co(num_orb_max,num_wann,num_atom,num_cell,ndim_kk_co_max)
! input: eigenvalues_co_rest(num_rest,nqdiv)
      subroutine cal_n_matrix
      
      use comlowh_mod

      implicit none
      include 'mpif.h'

      
      integer :: ik,irk,ibnd,jbnd,kbnd,iomega,
     $  ii,iio,icell,iatom,iorb,jorb,nspin,ispin
      double precision :: 
     $  n_electron_all
      complex*16 :: 
     $  g2mat(num_wann,num_wann),
     $  g3mat(num_wann,num_wann),
     $  gmat(num_wann,num_wann),
     $  g_lat_inv_co_temp(num_wann,num_wann),
     $  one_matrix_f(num_wann,num_wann),
     $  g_lat_inv_co(num_wann,num_wann,n_omega),
     $  mat_loc_ackw(num_orb_max,num_orb_max),
     $  tempmat1(num_orb_max,num_wann),
     $  exterp_g2(num_wann,num_wann),
     $  exterp_g3(num_wann,num_wann)


      ! one_matrix_f = 1
      one_matrix_f=0.0d0
      do ibnd=1, num_wann
        one_matrix_f(ibnd,ibnd)=1.0d0
      enddo

      ! decide nspin based on the value of iac.
      nspin=1
      do icell=1, num_cell
        do iatom=1, num_atom
          if (iac(iatom,icell) .lt. 0) then
            nspin=2
          endif
        enddo
      enddo

      n_matrix=0.0d0

       
      do ik=1, ndim_kk_co(me1)
        irk=n_mpi_kk_co(me1)+ik

        do ispin=1, nspin
          
          g_lat_inv_co=0.0d0

! g_lat_inv_co(omega) = frequency(iomega) - h_nonint_co(:,:,ik) + \sum_{atom} proj_renorm_co^\dagger * (sig(omega) - dc - I * gamma_selfenergy) * proj_renorm_co
!     embedding self-energy term: E(Sigma-DC)
          do iomega=1,n_omega
            
            do ibnd=1, num_wancell
              g_lat_inv_co(ibnd,ibnd,iomega)=frequency(iomega)
            enddo
            
            g_lat_inv_co(:,:,iomega)=
     $        g_lat_inv_co(:,:,iomega)-h_nonint_co(:,:,ik)          
            
            
            do icell=1, num_cell
              do iatom=1, num_atom

                if (ispin .eq. 1) then
                  ii= iac(iatom,icell)
                else
                  ii=-iac(iatom,icell)
                endif

! g_lat_inv_co(omega) = \sum_{atom} proj_renorm_co^\dagger * (sig(omega) - dc - I * gamma_selfenergy) * proj_renorm_co
                ! mat_loc_ackw = sig(omega) - dc - I * gamma_selfenergy
                mat_loc_ackw=0.0d0
                do iorb=1, num_orb(iatom)
                  do jorb=1, num_orb(iatom)

                    iio=ioac(iorb,jorb,ii)
                    if (iio .ne. 0) then
                      mat_loc_ackw(iorb,jorb)
     $                  =mat_loc_ackw(iorb,jorb)
     $                  +sig(iio,ii,iomega)
     $                  -dc(iio,ii)
                    endif
                    
                    if (iorb .eq. jorb) then
                      mat_loc_ackw(iorb,jorb)
     $                  =mat_loc_ackw(iorb,jorb)                    
     $                  -dcmplx(0.0d0, gamma_selfenergy) ! add broadening
                    endif ! (iorb .eq. jorb)

                  enddo ! jorb=1, num_orb(iatom)
                enddo ! iorb=1, num_orb(iatom)
                
                ! g_lat_inv_co+=proj_renorm_co^\dagger * mat_loc_ackw * proj_renorm_co
                tempmat1=0.0d0
                call zgemm('n','n',num_orb(iatom),num_wancell,
     $            num_orb(iatom),(1.d0,0.d0),
     $            mat_loc_ackw,num_orb_max,
     $            proj_renorm_co(1,1,iatom,icell,ik),num_orb_max,
     $            (0.d0,0.d0),
     $            tempmat1,num_orb_max)
                
                call zgemm('c','n',num_wancell,num_wancell,
     $            num_orb(iatom),(-1.d0,0.d0),
     $            proj_renorm_co(1,1,iatom,icell,ik),num_orb_max,
     $            tempmat1,num_orb_max,(1.d0,0.d0),
     $            g_lat_inv_co(1,1,iomega),num_wancell)

              enddo ! iatom=1, num_atom
            enddo ! icell=1, num_cell
          enddo !  iomega=1,n_omega
          
c$$$  nimp=0.0d0        
          
! Calculate the n_matrix:
! Use  gmat = (ef + g_lat_inv_co(:,:,n_omega))^{-1}  to calculate g2mat and g3mat, which will help with the tail at large momenta.
! Then use  gmat = (ef + g_lat_inv_co(:,:,iomega))^{-1} and g2mat, g3mat to calculate n_matrix.
! Finally add in n_matrix += \sum_{rest bands} 1/(1+exp(beta(eigenvalues_co_rest-ef))
          gmat=0.0d0        
          g_lat_inv_co_temp=g_lat_inv_co(:,:,n_omega)        
          
          do ibnd=1, num_wann
            g_lat_inv_co_temp(ibnd,ibnd)
     $        =g_lat_inv_co_temp(ibnd,ibnd)+ef
          enddo

c$$$  write(179,*), iomega, ik, irk

        ! gmat = (ef + g_lat_inv_co(:,:,n_omega)) ^{-1}
          call dcmplx_matinv(g_lat_inv_co_temp, gmat,
     $      num_wann,num_wann)
          
       ! g2mat = 0.5 * (gmat + gmat^\dagger) * frequency^2
          g2mat=(gmat+transpose(dconjg(gmat)))
     $      /2.0d0*frequency(n_omega)**2

       ! g3mat ?=? 0.5 * ((gmat - gmat^\dagger) -2/frequency) * frequency^3
          g3mat=(gmat-transpose(dconjg(gmat))
     $      -one_matrix_f*2.0/frequency(n_omega))
     $      /2.0d0*frequency(n_omega)**3
          
          exterp_g2=0.0d0
          exterp_g3=0.0d0        
          
          do iomega=1, n_omega
            
            gmat=0.0d0        
            g_lat_inv_co_temp=g_lat_inv_co(:,:,iomega)        
            
            do ibnd=1, num_wann
              g_lat_inv_co_temp(ibnd,ibnd)
     $          =g_lat_inv_co_temp(ibnd,ibnd)+ef
            enddo
c$$$  write(179,*), iomega, ik, irk

! gmat = (ef + g_lat_inv_co(:,:,iomega))^{-1}
            call dcmplx_matinv(g_lat_inv_co_temp, gmat,
     $        num_wann,num_wann)
            
! n_matrix = - beta^{-1} * ..... f(one_matrix_f, gmat, g2mat, g3mat)
            do ibnd=1,num_wann
              do jbnd=1,num_wann
                n_matrix(ibnd,jbnd,irk)
     $            =n_matrix(ibnd,jbnd,irk)
     $            -1.0d0/beta
     $            *(
     $            cdexp(-beta*frequency(iomega))
     $            *(
     $            gmat(ibnd,jbnd)
     $            -one_matrix_f(ibnd,jbnd)/frequency(iomega)
     $            -g2mat(ibnd,jbnd)/(frequency(iomega))**2
     $            -g3mat(ibnd,jbnd)/(frequency(iomega))**3                  
     $            )
     $            +
     $            cdexp(beta*frequency(iomega))
     $            *(
     $            dconjg(gmat(jbnd,ibnd))
     $            +one_matrix_f(ibnd,jbnd)/frequency(iomega)
     $            -g2mat(ibnd,jbnd)/(frequency(iomega))**2
     $            +g3mat(ibnd,jbnd)/(frequency(iomega))**3                  
     $            )
     $            )/nspin
                
                exterp_g2(ibnd,jbnd)
     $            =exterp_g2(ibnd,jbnd)
     $            -1.0d0/beta
     $            *(
     $            cdexp(-beta*frequency(iomega))
     $            *(
     $            -g2mat(ibnd,jbnd)/(frequency(iomega))**2
     $            )
     $            +
     $            cdexp(beta*frequency(iomega))
     $            *(
     $            -g2mat(ibnd,jbnd)/(frequency(iomega))**2
     $            )
     $            )/nspin
                
                exterp_g3(ibnd,jbnd)
     $            =exterp_g3(ibnd,jbnd)
     $            -1.0d0/beta
     $            *(
     $            cdexp(-beta*frequency(iomega))
     $            *(
     $            -g3mat(ibnd,jbnd)/(frequency(iomega))**3                  
     $            )
     $            +
     $            cdexp(beta*frequency(iomega))
     $            *(
     $            +g3mat(ibnd,jbnd)/(frequency(iomega))**3                  
     $            )
     $            )/nspin              
              enddo
            enddo
          enddo ! iomega=1, n_omega

! Now the loop over omega is done, but we have one more contribution to add to n_matrix, based on g2mat.
! n_matrix -= (-0.5 + 0.25 * g2mat * beta)/nspin
          do ibnd=1,num_wann
            do jbnd=1,num_wann
!     asymptotic interval
              n_matrix(ibnd,jbnd,irk)
     $          =n_matrix(ibnd,jbnd,irk)
     $          -1.0d0
     $          *(
     $          -one_matrix_f(ibnd,jbnd)/2.0d0            
     $          +g2mat(ibnd,jbnd)*beta/4.0d0
     $          )/nspin
              
              exterp_g2(ibnd,jbnd)
     $          =exterp_g2(ibnd,jbnd)
     $          -1.0d0
     $          *(
     $          +g2mat(ibnd,jbnd)*beta/4.0d0
     $          )/nspin            
            enddo ! jbnd=1,num_wann
          enddo !  ibnd=1,num_wann

        enddo ! ispin=1, nspin

c$$$        write(179,'(a, 2i5)')'g2mat',ik,irk
c$$$        do ibnd=1, num_wann
c$$$          write(179,'(10000(2f12.6,2x))')
c$$$     $      (g2mat(ibnd,jbnd), jbnd=1, num_wann)
c$$$        enddo
c$$$        write(179,'(a, 2i5)')'g3mat',ik,irk
c$$$        do ibnd=1, num_wann
c$$$          write(179,'(10000(2f12.6,2x))')
c$$$     $      (g3mat(ibnd,jbnd), jbnd=1, num_wann)
c$$$        enddo
c$$$        write(179,'(a, 2i5)')'exterp_g2',ik,irk
c$$$        do ibnd=1, num_wann
c$$$          write(179,'(10000(2f12.6,2x))')
c$$$     $      (exterp_g2(ibnd,jbnd), jbnd=1, num_wann)
c$$$        enddo
c$$$        write(179,'(a, 2i5)')'exterp_g3',ik,irk
c$$$        do ibnd=1, num_wann
c$$$          write(179,'(10000(2f12.6,2x))')
c$$$     $      (exterp_g3(ibnd,jbnd), jbnd=1, num_wann)
c$$$        enddo

! throw an error and crash if exterp_g3 is too big.
        if (maxval(cdabs(exterp_g3)) .gt. 1.0d-4) then
c$$$          write(179,'(2(f12.6, 2i5))') maxval(cdabs(exterp_g2)),
c$$$     $      maxloc(cdabs(exterp_g2)),
c$$$     $      maxval(cdabs(exterp_g3)),maxloc(cdabs(exterp_g3))                  
c$$$          write(179,*) 'higher order correction to the g is necessary'
          if (maswrk) then
            write(iun,'(2(f12.6, 2i5))') maxval(cdabs(exterp_g2)),
     $        maxloc(cdabs(exterp_g2)),
     $        maxval(cdabs(exterp_g3)),maxloc(cdabs(exterp_g3))        
            write(iun,*) 'higher order correction to the g is necessary'
          endif
          call ending
        endif

! n_matrix +=  fermi(beta(eigenvalues_co_rest-ef)), where eigenvalues_co_rest is a diagonal matrix and non-zero only in the P_rest sector.
        do ibnd=1, num_rest ! num_rest=num_bands-num_wann
          n_matrix(ibnd+num_wann,ibnd+num_wann,irk)
     $     =n_matrix(ibnd+num_wann,ibnd+num_wann,irk)
     $      +1.0d0
     $      /(dexp(beta*(eigenvalues_co_rest(ibnd,irk)-ef))+1.0d0)
        enddo

      enddo ! over momentum


      call mpi_allreduce_dcmplx(n_matrix,
     $  num_bands**2*nqdiv,mpi_sum,mpi_comm_world)
      

      end      

! cal_vw_matrix_co calculates vw_matrix_co and to_orig_bnd, which are used only by write_den_matrix. It also calculates eigenvalues_co_rest, which is used only by cal_n_matrix.
! vw_matrix_temp = [v_matrix * w_matrix_co, one_m_p] = [v_matrix * w_matrix_co, eigenvectors[1 - v_matrix^\dagger * v_matrix]]
! vw_matrix_co = v_matrix * w_matrix_co + vw_matrix_temp * eigenvectors[P_rest [vw_matrix_temp^\dagger * eigenvalues * vw_matrix_temp] P_rest]
! P_rest projects out the num_rest=num_bands-num_wann sector
! to_orig_bnd is determined by finding the maximum element within specific vectors in vw_matrix_co.
! eigenvalues_co_rest is the eigenvalues  of P_rest [vw_matrix_temp^\dagger * eigenvalues * vw_matrix_temp] P_rest.
! input: v_matrix
! input: eigenvalues
! input: w_matrix_co
! ----------------
! output: eigenvalues_co_rest = eigenvalues[ P_r [vw_matrix_temp^\dagger * eigenvalues * vw_matrix_temp] P_r ]
! output: vw_matrix_co = v_matrix(q)  P_w w_matrix_co(q) P_w + eigenvectors[1 - v_matrix^\dagger P_w v_matrix] P_r * eigenvectors [ P_r vw_matrix_temp^\dagger * eigenvalues * vw_matrix_temp P_r ]
! output: to_orig_bnd.  jbnd=to_orig_bnd(ibnd,k) is an integer specifying which vector in vw_matrix_co(:,jbnd,k) matches most strongly to the basis element (band) specified by ibnd.
! vw_matrix_temp = v_matrix(q)  P_w w_matrix_co(q) P_w + eigenvectors[1 - v_matrix^\dagger P_w v_matrix] P_r
! v_matrix is from wannier.dat, and w_matrix_co is calculated by h_nonint_co.
      subroutine cal_vw_matrix_co

      use comlowh_mod

      implicit none
      include 'mpif.h'
      
      integer :: ik,irk,tempind(1),bndind_rest(num_rest),ii,
     $  ibnd,jbnd,kbnd,tempvec_int(num_bands),
     $  tempint(1)
! num_rest=num_bands-num_wann

      double precision :: tempvec_dble(num_bands)
      complex*16 :: vw_matrix_temp(num_bands,num_bands),
     $  hmat_rest(num_rest,num_rest),tempvec(num_bands),
     $  one_m_p(num_bands,num_bands)
      
      do ik=1,ndim_kk_co(me1)

        vw_matrix_temp=0.0d0
        irk=n_mpi_kk_co(me1)+ik
        
! vw_matrix_co(q) = v_matrix(q)  P_w w_matrix_co(q) P_w
! P_w projects the wannier subspace out of the bands.
! v_matrix is from wannier.dat, and w_matrix_co is calculated by h_nonint_co.
! vw_matrix_co is num_bands x num_bands, v_matrix is  num_bands x  num_wann, and w_matrix_co is num_wann x num_wann
        call zgemm('n','n',num_bands,num_wann,
     $    num_wann,(1.d0,0.d0),
     $    v_matrix(1,1,irk),num_bands,
     $    w_matrix_co(1,1,ik),num_wann,(0.d0,0.d0),
     $    vw_matrix_co(1,1,irk),num_bands)
                
! one_m_p = 1
        one_m_p=0.0d0
        do ibnd=1, num_bands
          one_m_p(ibnd,ibnd)=1.0d0
        enddo

! one_mp = 1 - v_matrix P_w v_matrix^\dagger
        ! projection operator to the wannier space
        call zgemm('n','c',num_bands,num_bands,
     $    num_wann,(-1.d0,0.d0),
     $    v_matrix(1,1,irk),num_bands,
     $    v_matrix(1,1,irk),num_bands,(1.d0,0.d0),
     $    one_m_p,num_bands)

! tempvec_dble and one_m_p are the eigenvalues and eigenvectors of 1 - v_matrix P_w v_matrix^\dagger
! one_m_p = eigenvectors[1 - v_matrix P_w v_matrix^\dagger]
! tempvec_dble = eigenvalues[1 - v_matrix P_w v_matrix^\dagger]
        call hermitianeigen_cmplxdouble
     $    (num_bands,tempvec_dble,one_m_p)
c$$$        write(179,'(a, i5, 10000f12.6)')
c$$$     $    'one_m_p eigenvalues', irk, tempvec_dble(1:num_wann)
c$$$        write(179,'(a, i5, 10000f12.6)')
c$$$     $    'one_m_p eigenvalues', irk,
c$$$     $    tempvec_dble((num_wann+1):num_bands)


        vw_matrix_temp=vw_matrix_co(:,:,irk)
! vw_matrix_temp = v_matrix(q)  P_w w_matrix_co(q) P_w

        do ibnd=1, num_rest ! num_rest=num_bands-num_wann
          vw_matrix_temp(:,ibnd+num_wann)=one_m_p(:,ibnd+num_wann)
        enddo
! vw_matrix_temp = v_matrix(q)  P_w w_matrix_co(q) P_w + eigenvectors[1 - v_matrix^\dagger P_w v_matrix] P_r
! P_r = 1 - P_w

! is_orthonormal_dcmplx calculates vw_matrix_temp^\dagger vw_matrix_temp and checks that it is equal to one. If it is not, terminate execution.
        call is_orthonormal_dcmplx
     $    (vw_matrix_temp,num_bands,
     $    vw_matrix_temp,num_bands,
     $    num_bands,num_bands,num_bands,1,
     $    'basis_after_orthonormal_basis_construction')
        
! hmat_rest = P_r [vw_matrix_temp^\dagger * eigenvalues * vw_matrix_temp] P_r, where P_r = 1 - P_w projects out the num_rest=num_bands-num_wann sector
!     disentangle the rest bands
        hmat_rest=0.0d0
        do ibnd=1, num_rest ! num_rest=num_bands-num_wann
          do jbnd=1, num_rest
            do kbnd=1, num_bands
              hmat_rest(ibnd,jbnd)=hmat_rest(ibnd,jbnd)
     $          +dconjg(vw_matrix_temp
     $          (kbnd,ibnd+num_wann))
     $          *eigenvalues(kbnd,irk)
     $          *vw_matrix_temp
     $          (kbnd,jbnd+num_wann)
            enddo
          enddo
        enddo
        
! eigenvalues_co_rest = eigenvalues[ P_r [vw_matrix_temp^\dagger * eigenvalues * vw_matrix_temp] P_r ]
! hmat_rest = eigenvectors[ P_r [vw_matrix_temp^\dagger * eigenvalues * vw_matrix_temp] P_r ]
        call hermitianeigen_cmplxdouble(num_rest,
     $    eigenvalues_co_rest(:,irk),hmat_rest)
        

! vw_matrix_co = v_matrix(q)  P_w w_matrix_co(q) P_w + vw_matrix_temp * hmat_rest
! vw_matrix_co = v_matrix(q)  P_w w_matrix_co(q) P_w + eigenvectors[1 - v_matrix^\dagger P_w v_matrix] P_r * eigenvectors [ P_r vw_matrix_temp^\dagger * eigenvalues * vw_matrix_temp P_r ]
! vw_matrix_temp = v_matrix(q)  P_w w_matrix_co(q) P_w + eigenvectors[1 - v_matrix^\dagger P_w v_matrix] P_r
! hmat_rest = eigenvectors [ P_r vw_matrix_temp^\dagger * eigenvalues * vw_matrix_temp P_r ]
        call zgemm('n','n',num_bands,num_rest,num_rest,
     $    (1.d0,0.d0),vw_matrix_temp(1,num_wann+1),num_bands,
     &    hmat_rest,num_rest,(1.d0,0.d0),
     $    vw_matrix_co(1,num_wann+1,irk),num_bands)
        
! check that vw_matrix_co^\dagger vw_matrix_co = 1, and if not terminate execution.
        call is_orthonormal_dcmplx
     $    (vw_matrix_co(1,1,irk),num_bands,
     $    vw_matrix_co(1,1,irk),num_bands,
     $    num_bands,num_bands,num_bands,1,
     $    'final vw_matrix')

      enddo ! ik=1,ndim_kk_co(me1)

      call mpi_allreduce_dcmplx(vw_matrix_co,
     $  num_bands**2*nqdiv,mpi_sum,mpi_comm_world)
      call mpi_allreduce_dble(eigenvalues_co_rest,
     $  num_rest*nqdiv,mpi_sum,mpi_comm_world)

! check again that vw_matrix_co is unitary, and if not terminate execution.
      if (maswrk) then
        do irk=1, nqdiv
          call is_orthonormal_dcmplx
     $      (vw_matrix_co(1,1,irk),num_bands,
     $      vw_matrix_co(1,1,irk),num_bands,
     $      num_bands,num_bands,num_bands,1,
     $      'ff vw_matrix')
        enddo
      endif

! The whole point of this next block is to calculate to_orig_bnd, which is used ONLY in write_den_matrix.
! jbnd=to_orig_bnd(ibnd,k) is an integer specifying which vector in vw_matrix_co(:,jbnd,k) matches most strongly to the basis element (band) specified by ibnd.
! this looks at the vector vw_matrix_co(:,ibnd,k).  Based on the which element of this vector has the largest magnitude, it sets to_orig_bnd(ibnd,k) equal to the index of that element.  The bits about tempvec_int guarantee that orig_bnd does not have duplicates.
      to_orig_bnd=0
      
      do irk=1, nqdiv
        tempvec_int=0
        
        do ibnd=1, num_bands
          tempvec_dble=0.0d0

          do jbnd=1, num_bands

            tempvec_dble(jbnd)
     $        =cdabs(vw_matrix_co(jbnd,ibnd,irk))
            if (tempvec_int(jbnd) .eq. 1) then
              tempvec_dble(jbnd)=0.0d0
            endif

          enddo ! jbnd=1, num_bands
          
          tempint=maxloc(tempvec_dble)
          to_orig_bnd(ibnd,irk)=tempint(1)
          tempvec_int(tempint(1))=1

        enddo ! ibnd=1, num_bands

      enddo !  irk=1, nqdiv

c$$$      if (maswrk) then
c$$$        do irk=1, nqdiv      
c$$$          do ibnd=1, num_bands
c$$$            jbnd=to_orig_bnd(ibnd,irk)
c$$$            if (ibnd .le. num_wann) then
c$$$              write(iun,'(3i5, 3f12.6)'), irk,
c$$$     $          jbnd,ibnd,eigenvalues(jbnd,irk),
c$$$     $          eigenvalues_co_all(ibnd,irk),
c$$$     $          eigenvalues(jbnd,irk)
c$$$     $          -eigenvalues_co_all(ibnd,irk)
c$$$            else
c$$$              write(iun,'(3i5, 3f12.6)'), irk,
c$$$     $          jbnd,ibnd,eigenvalues(jbnd,irk),
c$$$     $          eigenvalues_co_rest(ibnd-num_wann,irk),
c$$$     $          eigenvalues(jbnd,irk)
c$$$     $          -eigenvalues_co_rest(ibnd-num_wann,irk)              
c$$$     $          
c$$$            endif
c$$$          enddo
c$$$        enddo
c$$$      endif
      end
      
      subroutine cal_max_or_min_ind(inarray, inarray_n,
     $  outind, outind_n, flag)
      implicit none
      integer, intent(in) :: inarray_n, outind_n
      double precision, intent(in) :: inarray(inarray_n)
      character*3, intent(in) :: flag      
      integer, intent(out) :: outind(outind_n)

      integer :: tempind(1),ii      
      double precision :: inarray_copy(inarray_n), bound

      inarray_copy=inarray
      outind=0

      if (flag .eq. 'min') then

        bound=maxval(inarray_copy)+1.0d0
        do ii=1, outind_n
          tempind=minloc(inarray_copy)
          inarray_copy(tempind(1))=bound
          outind(ii)=tempind(1)
        enddo
      elseif (flag .eq. 'max') then

        bound=minval(inarray_copy)-1.0d0
        do ii=1, outind_n
          tempind=maxloc(inarray_copy)
          inarray_copy(tempind(1))=bound
          outind(ii)=tempind(1)
        enddo
      endif
      end

      subroutine find_bnd_in_win(eigenval,nk,nbnd,emin,emax,
     $  den_mode,is_spinorbit,nbnd_k,ind_bnd_k)
! find_bnd_in_win creates ind_bnd_k which says which bands are inside the energy window [emin,emax] for at least one k-point, and nbnd_k which is the count of bands inside the window.  Even though ind_bnd_k and nbnd_k are indexed over k-points, they have the same value for all k-points.
! eigenval <=> eigenvalues_fi_all
! nk <=> numk
! nbnd <=> num_wancell = num_wann
! emin <=> ewin_min
! emax <=> ewin_max
! den_mode = 'bnd'
! output: nbnd_k <=> nbnd_in_win
! output: ind_bnd_k <=> bnd_in_win
      implicit none
      include 'mpif.h'
      integer, intent(in) :: nk, nbnd, is_spinorbit
      double precision, intent(in) :: eigenval(nbnd,nk),emin,emax
      character*3, intent(in)  :: den_mode

      integer, intent(out) :: nbnd_k(nk)
      integer, intent(out) :: ind_bnd_k(nbnd,nk)

      integer :: irk,ibnd,ik(1),ibnd_min,ibnd_max
      
! If emin < eigenval(ibnd,irk) < emax then we mark this band as being inside the window [in ind_bnd_k], and we increment the count of bands inside the window [nbnd_k]. These are temporary values; the final values of ind_bnd_k and nbnd_k will be decided at the end.
      nbnd_k=0
      ind_bnd_k=0
      do irk=1, nk
        do ibnd=1, nbnd
          if ((eigenval(ibnd,irk) .gt. emin) .and.
     $      (eigenval(ibnd,irk) .lt. emax)) then
            nbnd_k(irk)=nbnd_k(irk)+1
            ind_bnd_k(nbnd_k(irk),irk)=ibnd
          endif
        enddo
      enddo

! den_mode='bnd'
      if (den_mode.eq. 'bnd') then

! ik is the index of k-point which has the most bands inside the window.
        ik=maxloc(nbnd_k)
        ibnd_min=ind_bnd_k(1,ik(1))
        ibnd_max=ind_bnd_k(nbnd_k(ik(1)),ik(1))
c$$$        write(179,*) 'start', ik(1),ibnd_min,ibnd_max

! This loop discovers minimum and  maximum band indices which are inside the window for at least one k-point: ibnd_min, ibnd_max.
        do irk=1, nk
          do ibnd=1, ibnd_min-1
            if (dabs(eigenval(ibnd,irk)-eigenval(ibnd_min,irk))
     $        .lt. 1.0d-4) then
              ibnd_min=ibnd
            endif
          enddo

          do ibnd=ibnd_max+1, nbnd
            if (dabs(eigenval(ibnd,irk)-eigenval(ibnd_max,irk))
     $        .lt. 1.0d-4) then
c$$$              write(179,'(3i5,2f12.6)') irk,ibnd,ibnd_max,
c$$$     $          eigenval(ibnd,irk),eigenval(ibnd_max,irk)
              ibnd_max=ibnd              
            endif
          enddo
        enddo ! irk=1, nk

! Adjust ibnd_min, ibnd_max for spin orbit.
        if (is_spinorbit==1) then
          if (mod(ibnd_min,2)==0) then
            ibnd_min=ibnd_min-1
          endif
          if (mod(ibnd_max,2)==1) then
            ibnd_max=ibnd_max+1
          endif
        endif

        ind_bnd_k=0
        nbnd_k=0
        do irk=1, nk
          nbnd_k(irk)=ibnd_max-ibnd_min+1
          do ibnd=ibnd_min, ibnd_max
            ind_bnd_k(ibnd-ibnd_min+1,irk)=ibnd
          enddo
        enddo

      endif ! (den_mode.eq. 'bnd')

      end ! subroutine find_bnd_in_win
      

! write_den_matrix writes out wannier_den_matrix.dat, which is used in rspflapw.
! wannier_den_matrix.dat can be in either binary or hdf5 format.
! It writes:
!       - "temperature"
!       - "nqdiv" - the number of k-points
!       - "n_dmft_bnd" from nbnd_in_win_co [called 'num_band_in' when saved from comrisb, and read in by rspflapw as 'n_dmft_bnd'] - should be the number of bands
!       - "n_dmft_wan" from nbnd_in_win_co 
!       - "include_band_new"
!               include_band_new comes from nbnd_in_win_co, bnd_in_win_co, to_orig_bnd, and include_bands.
!       - "vw_matrix_new"
!               vwmatrix_new is a subsector of vw_matrix_co, with the sector specified by nbnd_in_win_co, to_orig_bnd, bnd_in_win_co. vwmatrix_new is renormalized so that it satisfies vwmatrix_new^\dagger vwmatrix_new = 1.
!       - 'n_matrix_new' - the density matrix
!                n_matrix_new equal is a subsector of n_matrix, with the sector specified by nbnd_in_win_co, bnd_in_win_co.
      subroutine write_den_matrix

      use comlowh_mod

#ifdef USE_HDF5
      use hdf5io_mod
      use hdf5
#endif


      implicit none
      include 'mpif.h'
#ifdef USE_HDF5
      integer :: ng(7)
      integer :: nv(7)
      integer(hid_t) :: fid_dat = -1
      integer :: h5err
#endif


      integer :: 
     $  nbnd_in_win_co_max,
     $  irk,ibnd,jbnd,kbnd,abnd,bbnd,ik
      double precision :: n_electron_all

      integer, allocatable :: include_band_new(:,:)
      complex*16, allocatable :: n_matrix_new(:,:,:),
     $  vw_matrix_new(:,:,:),proj_in(:,:),proj_out(:,:)

!  nbnd_in_win_co=nbnd_in_win(1), bnd_in_win_co=bnd_in_win(:,1)

      ! This decides include_band_new, based on nbnd_in_win_co, bnd_in_win_co, to_orig_bnd, and include_bands. include_band_new will be saved to wannier_den_matrix.dat
      nbnd_in_win_co_max=maxval(nbnd_in_win_co)
      allocate(include_band_new(nbnd_in_win_co_max,nqdiv))
      include_band_new=0
      do irk=1, nqdiv
        do ibnd=1, nbnd_in_win_co(irk)
          jbnd=bnd_in_win_co(ibnd,irk)
          kbnd=to_orig_bnd(jbnd,irk)
          include_band_new(ibnd,irk)=include_bands(kbnd)
! include_bands is from wannier.dat
        enddo
      enddo      
      
      allocate(n_matrix_new(nbnd_in_win_co_max,nbnd_in_win_co_max,
     $  nqdiv))
      n_matrix_new=0.0d0
      
      allocate(vw_matrix_new(nbnd_in_win_co_max,nbnd_in_win_co_max,
     $  nqdiv))
      vw_matrix_new=0.0d0        
      
      do ik=1, ndim_kk_co(me1)
        irk=n_mpi_kk_co(me1)+ik        
        
        allocate(proj_in(nbnd_in_win_co(irk),nbnd_in_win_co(irk)))
        proj_in=0.0d0
        allocate(proj_out(nbnd_in_win_co(irk),nbnd_in_win_co(irk)))
        proj_out=0.0d0          

        ! This sets proj_in equal to a subsector of vw_matrix_co, with the sector specified by nbnd_in_win_co, to_orig_bnd, bnd_in_win_co.
% proj_in = U_{tob} P_E  v_{wmc}(q,i,j) P_E
        do ibnd=1, nbnd_in_win_co(irk)
          kbnd=bnd_in_win_co(ibnd,irk)
          abnd=to_orig_bnd(kbnd,irk)            
          do jbnd=1, nbnd_in_win_co(irk)
            bbnd=bnd_in_win_co(jbnd,irk)
            
            proj_in(ibnd,jbnd)
     $        =proj_in(ibnd,jbnd)
     $        +vw_matrix_co(abnd,bbnd,irk)
          enddo
        enddo
        
% v_{wmn} = proj_in (proj_in proj_in^\dagger)^{-1}
       ! projector_renormalization returns proj_out = proj_in / \sqrt{proj_in^\dagger proj_in}.  proj_out has the property that proj_out^\dagger proj_out = 1.
        call projector_renormalization
     $    (nbnd_in_win_co(irk),nbnd_in_win_co(irk),proj_in,proj_out)
       
        ! This sets vw_matrix_new = proj_out. vw_matrix_new is saved in wannier_den_matrix.
        do ibnd=1, nbnd_in_win_co(irk)
          do jbnd=1, nbnd_in_win_co(irk)
            vw_matrix_new(ibnd,jbnd,irk)
     $        =vw_matrix_new(ibnd,jbnd,irk)
     $        +proj_out(ibnd,jbnd)
          enddo
        enddo

        deallocate(proj_in)
        deallocate(proj_out)
        
! n_matrix_new = P_E n_matrix(q,i,j)} P_E
! This sets n_matrix_new equal to a subsector of n_matrix, with the sector specified by nbnd_in_win_co, bnd_in_win_co. n_matrix_new is saved in wannier_den_matrix.dat.
        do ibnd=1, nbnd_in_win_co(irk)
          abnd=bnd_in_win_co(ibnd,irk)            
          do jbnd=1, nbnd_in_win_co(irk)
            bbnd=bnd_in_win_co(jbnd,irk)              
            n_matrix_new(ibnd,jbnd,irk)
     $        =n_matrix_new(ibnd,jbnd,irk)
     $        +n_matrix(abnd,bbnd,irk)
          enddo
        enddo        
      enddo
      
      call mpi_allreduce_dcmplx(n_matrix_new,
     $  nbnd_in_win_co_max**2*nqdiv,mpi_sum,mpi_comm_world)      
      
      call mpi_allreduce_dcmplx(vw_matrix_new,
     $  nbnd_in_win_co_max*nbnd_in_win_co_max*nqdiv,
     $  mpi_sum,mpi_comm_world)

#ifdef USE_HDF5
      nv=0
      call fopen('wannier_den_matrix.dat',fid_dat)
      call fwrite(MPI_COMM_WORLD,temperature,"temperature",fid_dat)
      call fwrite(MPI_COMM_WORLD,nqdiv,"nqdiv",fid_dat)
      ng(1)=size(nbnd_in_win_co,1)
      call fwrite(MPI_COMM_WORLD,nbnd_in_win_co,ng(1),nv(1),
     $                "n_dmft_bnd",fid_dat)
      ng(1)=size(nbnd_in_win_co,1)
      call fwrite(MPI_COMM_WORLD,nbnd_in_win_co,ng(1),nv(1),
     $                "n_dmft_wan",fid_dat)     
      ng(1)=size(include_band_new,1)
      ng(2)=size(include_band_new,2)
      call fwrite(MPI_COMM_WORLD,include_band_new,ng,nv,
     $                "include_band_new",fid_dat)
      ng(1)=size(vw_matrix_new,1)
      ng(2)=size(vw_matrix_new,2)
      ng(3)=size(vw_matrix_new,3)
      call fwrite(MPI_COMM_WORLD,vw_matrix_new,ng,nv,
     $  "vw_matrix_new",fid_dat)
      if (maswrk) write(iun,*)'vw_matrix_new',vw_matrix_new
      if (maswrk) write(iun,*)'size',ng(1),ng(2),ng(3)
      ng(1)=size(n_matrix_new,1)
      ng(2)=size(n_matrix_new,2)
      ng(3)=size(n_matrix_new,3)
      call fwrite(MPI_COMM_WORLD,n_matrix_new,ng,nv,
     $ "n_matrix_new",fid_dat)
      if (maswrk) write(iun,*)'n_matrix_new',n_matrix_new
      if (maswrk) write(iun,*)'size',ng(1),ng(2),ng(3)

      call fclose(fid_dat)
#else
      if (maswrk) then
        open(unit=305, file='wannier_den_matrix.dat',
     $    form='unformatted')
        write(305) temperature
        write(305) nqdiv
        write(305) nbnd_in_win_co
        write(305) nbnd_in_win_co
        write(305) include_band_new
        write(305) vw_matrix_new
        write(305) n_matrix_new
        close(305)
        
        write(iun,*) 'num_band information'
        write(iun,'(a, 10000i5)') 'bands in the window',
     $    nbnd_in_win_co_max, nbnd_in_win_co
        do irk=1, nqdiv
          do ibnd=1, nbnd_in_win_co(irk)
            jbnd=bnd_in_win_co(ibnd,irk)
            kbnd=to_orig_bnd(jbnd,irk)            
            write(iun,'(4i5, 2f12.6)')
     $        ibnd, jbnd,kbnd,include_band_new(ibnd,irk),
     $        eigenvalues_co_all(jbnd,irk), eigenvalues(kbnd,irk)
          enddo
        enddo
      endif
#endif

      deallocate(include_band_new)
      deallocate(n_matrix_new)
      deallocate(vw_matrix_new)                
      
      end      
