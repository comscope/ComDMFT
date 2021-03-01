      subroutine read_ini
      use manager_mod
      use wannier90_mod
      use parallel_mod
      use solid_mod
      use atom_mod
      use units_mod
      use vertex_mod
#ifdef USE_HDF5
      use hdf5io_mod
      use hdf5
#endif

      implicit none
      include 'mpif.h'
      integer :: ii,jj,i,
     $  cnt1,cnt2,iwan
      double precision :: real_lattice_temp(3,3),recip_lattice_temp(3,3)
#ifdef USE_HDF5
      integer(hid_t) :: fid_gwl = -1
      integer :: h5err
      integer :: ng(4)
      integer :: nv(4)
#endif
#ifdef USE_HDF5

      if(maswrk) then
        open(1,file='comcoulomb.ini',status='old')
        read(1,'(a)') rstpath
        read(1,'(a)') wanpath
        read(1,*) nproc_tau
        read(1,*) nproc_k
        close(1)
      endif
      if(goparr) then
        call mpi_bcast(rstpath,1000,
     $    mpi_character,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(wanpath,1000,
     $    mpi_character,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nproc_tau,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nproc_k,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      endif
      call finit()
      call fopen(adjustl(trim(wanpath))//'/wannier.dat',fid_gwl)

      nv=0

      ng(1)=size(real_lattice_temp,1)
      ng(2)=size(real_lattice_temp,2)
      call fread(MPI_COMM_WORLD,real_lattice_temp,ng,nv,"rbas",fid_gwl)
      ng(1)=size(recip_lattice_temp,1)
      ng(2)=size(recip_lattice_temp,2)
      call fread(MPI_COMM_WORLD,recip_lattice_temp,ng,nv,"gbas",fid_gwl)
      call fread(MPI_COMM_WORLD,num_bands,"num_bands",fid_gwl)
      call fread(MPI_COMM_WORLD,num_wann,"num_wann",fid_gwl)
      call fclose(fid_gwl)
#endif
      if(maswrk) then
        open(1,file='comcoulomb.ini',status='old')
        read(1,'(a)') rstpath
        read(1,'(a)') wanpath
        read(1,*) nproc_tau
        read(1,*) nproc_k
#ifdef USE_HDF5

#else
        open(unit=10, file=adjustl(trim(wanpath))//'/wannier.dat',
     $    form='unformatted')
        read(10) real_lattice_temp(3,3)
        read(10) recip_lattice_temp(3,3)      
        read(10) num_bands      !!
        read(10) num_wann       !!
        close(10)
#endif        
        read(1,*)projemin
        read(1,*)projemax
        projemax=projemax*2.d0/evolt
        projemin=projemin*2.d0/evolt
        read(1,*)include_gamma_w
        read(1,*)u_mode                
        if ((u_mode .eq. 'wnds') .or. (u_mode .eq. 'gtop') .or.
     $    (u_mode .eq. 'edmf') .or. (u_mode .eq. 'bnse')) then
          read(1,*) num_orb_cor
          if(allocated(wan_correlated))deallocate(wan_correlated)
          allocate(wan_correlated(num_orb_cor))
          read(1,*) (wan_correlated(iwan), iwan=1, num_orb_cor)
          do ii=1, num_orb_cor
            if (wan_correlated(ii) .gt. num_wann) then
              write(iun,*) 'wan_correlated wrong'
              call flush(iun)
              call ending
            endif
          enddo
          enpj_emin=-80.0d0
          enpj_emax=80.0d0

        elseif (u_mode .eq. 'enpj') then
          read(1,*)enpj_emin
          read(1,*)enpj_emax          
          enpj_emin=enpj_emin*2.d0/evolt
          enpj_emax=enpj_emax*2.d0/evolt
          num_orb_cor=num_wann
          if(allocated(wan_correlated))deallocate(wan_correlated)
          allocate(wan_correlated(num_orb_cor))
          do iwan=1, num_orb_cor
            wan_correlated(iwan)=iwan
          enddo
        endif

        if ((u_mode .eq. 'wnds') .or. (u_mode .eq. 'gtop') .or.
     $    (u_mode .eq. 'enpj').or. (u_mode .eq. 'bnse')) then        
          read(1,*)p_low_mode
          read(1,*)p_local_k
          read(1,*)dist_rc        
          read(1,*)u_linear
        elseif (u_mode .eq. 'edmf') then
          p_low_mode=1
          p_local_k=.false.
          dist_rc=30.0d0
          u_linear=.false.
        endif

        ewald=2.0d0
        ewald0=ewald        
        r_ref=0.0d0
        p_local_r=.false.


        num_orb_iter=num_wann-num_orb_cor
        if (num_orb_iter .gt. 0) then
          if(allocated(wan_itinerant))
     $      deallocate(wan_itinerant)
          allocate(wan_itinerant(num_orb_iter))
          cnt2=0
          do ii=1, num_wann
            cnt1=0
            do jj=1, num_orb_cor
              if (ii .eq. wan_correlated(jj)) then
                cnt1=cnt1+1
              endif
            enddo
            if (cnt1 .eq. 0) then
              cnt2=cnt2+1
              wan_itinerant(cnt2)=ii
            endif
          enddo
        endif

        num_orb_proj=num_orb_cor+num_orb_iter
        
        close(1)

c$$$        write(iun,*) adjustl(trim(rstpath))
c$$$        write(iun,*) adjustl(trim(wanpath))
c$$$        write(iun,*) nproc_tau
c$$$        write(iun,*) nproc_k        
c$$$
c$$$        write(iun,*) num_bands      !!
c$$$        write(iun,*) num_wann       !!
c$$$        
c$$$        write(iun,*)projemin*rydberg
c$$$        write(iun,*)projemax*rydberg
c$$$        write(iun,*)include_gamma_w
c$$$        write(iun,*)u_mode                
c$$$        write(iun,*) num_orb_cor
c$$$        write(iun,*) (wan_correlated(iwan), iwan=1, num_orb_cor)
c$$$        write(iun,*)enpj_emin
c$$$        write(iun,*)enpj_emax          
c$$$        write(iun,*)p_low_mode
c$$$        write(iun,*)p_local_k
c$$$        write(iun,*)dist_rc        
c$$$        write(iun,*)u_linear

        ewald=2.0d0
        ewald0=ewald        
        r_ref=0.0d0
        p_local_r=.false.        
      endif

      if(goparr) then
        call mpi_bcast(rstpath,1000,
     $    mpi_character,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(wanpath,1000,
     $    mpi_character,0,MPI_COMM_WORLD,mpierr)                
        call mpi_bcast(nproc_tau,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nproc_k,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(enpj_emin,1,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(enpj_emax,1,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(u_mode,4,mpi_character,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(p_low_mode,1,mpi_integer,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(p_local_r,1,mpi_logical,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(p_local_k,1,mpi_logical,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(u_linear,1,mpi_logical,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast
     $    (projemin,1,mpi_double_precision,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast
     $    (projemax,1,mpi_double_precision,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(ewald,1,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(ewald0,1,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(dist_rc,1,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(include_gamma_w,1,mpi_logical,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(r_ref,3,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(num_orb_cor,1,mpi_integer,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(num_orb_iter,1,mpi_integer,
     $    0,MPI_COMM_WORLD,mpierr)
        num_orb_proj=num_orb_cor+num_orb_iter
        call mpi_bcast(num_orb_proj,1,mpi_integer,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(num_wann,1,mpi_integer,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(num_bands,1,mpi_integer,
     $    0,MPI_COMM_WORLD,mpierr)
        if (.not. maswrk)then
          if(allocated(wan_correlated))
     $      deallocate(wan_correlated)
          allocate(wan_correlated(num_orb_cor))
        endif
        call mpi_bcast(wan_correlated,num_orb_cor,mpi_integer,
     $    0,MPI_COMM_WORLD,mpierr)
        if (num_orb_iter .gt. 0)  then
          if (.not. maswrk)then
            if(allocated(wan_itinerant)) 
     $        deallocate(wan_itinerant)
            allocate(wan_itinerant(num_orb_iter))
          endif
          call mpi_bcast(wan_itinerant,num_orb_iter,mpi_integer,
     $      0,MPI_COMM_WORLD,mpierr)
        endif
      endif       

      end
