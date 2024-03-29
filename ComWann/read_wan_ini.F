      subroutine read_comwann_ini
      use comwann_mod
      implicit none
      include 'mpif.h'

      character*300 :: ubi_txt
      
      if(maswrk) then
        open(11,file='comwann.ini',status='old')
        read(11,'(a)') rstpath
        read(11,'(a)') ubi_txt
        if (adjustl(trim(ubi_txt)) .eq. 'dft') then
          ubi_num=1
          ubi_0='dft'
c$$$  rmode=1
c$$$  mt_fac=1.0d0
        elseif (adjustl(trim(ubi_txt)) .eq. 'qp') then
          ubi_num=4
          ubi_0='qp'
c$$$  rmode=2
c$$$  mt_fac=0.999d0
        endif
        write(993, *) 'ubi_num=', ubi_num
        write(993,*) 'ubi_0=',ubi_0        

        read(11,*) dis_win_max
        read(11,*) dis_win_min
        read(11,*) dis_froz_max
        read(11,*) dis_froz_min
        read(11,*) num_iter
        read(11,*) dis_num_iter        
        read(11,*) writewan
        read(11,*) cut_low
        read(11,*) cut_froz
        read(11,*) cut_total
        read(11,*) rmode
        read(11,*) mt_fac
        read(11,*) dis_froz_max_fac
        write(993,*) 'dis_win_max=',dis_win_max
        write(993,*) 'dis_win_min=',dis_win_min
        write(993,*) 'dis_froz_max=',dis_froz_max
        write(993,*) 'dis_froz_min=',dis_froz_min
        write(993,*) 'num_iter=',num_iter
        write(993,*) 'dis_num_iter=',dis_num_iter        
        write(993, *) 'write_wan=', writewan
        
c$$$  cut_low=0.4d0
c$$$  cut_froz=0.15d0
c$$$  cut_total=0.0d0
        write(993, *) 'cut_low=', cut_low
        write(993, *) 'cut_froz=', cut_froz
        write(993, *) 'cut_total=', cut_total
        write(993, *) 'rmode=', rmode
        write(993, *) 'mt_fac=', mt_fac
        write(993, *) 'dis_froz_max_fac=', dis_froz_max_fac
      endif
      if(goparr) then
        call mpi_bcast(rstpath,1000,
     $    mpi_character,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(ubi_num,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(ubi_0,3,mpi_character,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(dis_froz_min,1, mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(dis_froz_max,1, mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(dis_win_min,1, mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)        
        call mpi_bcast(dis_win_max,1, mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(num_iter,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(dis_num_iter,1,mpi_integer,
     $    0,MPI_COMM_WORLD,mpierr)                        
        call mpi_bcast(writewan,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(cut_low,1,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(cut_froz,1,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)                
        call mpi_bcast(cut_total,1,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(rmode,1,mpi_integer,0,
     $    MPI_COMM_WORLD,mpierr)                                        
        call mpi_bcast(mt_fac,1,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(dis_froz_max_fac,1,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)                                        
      endif
      if (.not. ((ubi_num .eq. 1) .or. (ubi_num .eq. 4))) then
        write(993,*) 'Please enter correct cal_mode either dft or qp'
        call ending_qft
      endif

      close(11)
      end


      subroutine write_wannier_win(ispin)
      use comwann_mod

      implicit none
      integer, intent(in) :: ispin
      integer :: iatom,ii,nlines,io
      logical :: file_exists
!     double precision :: eval1,eval2,ee

      character*300 :: tempchar
      character*500 :: tempstr

      if (maswrk) then


        inquire(file='kpath',exist=file_exists)
        if (file_exists) then
          nlines = 0 
          OPEN (1, file = 'kpath')
          DO
            READ(1,*,iostat=io)
            IF (io/=0) EXIT
            nlines = nlines + 1
          END DO
          CLOSE (1)
        endif
        
        open(unit=992, file=adjustl(trim(seed_name))//'.win')
        write(992,'(a)') "begin Projections"
        write(992,'(a)') "random"
        write(992,'(a)') "end Projections"
        
        write(992,'(a)') "begin atoms_cart"
        write(992,'(a)') "Ang"
        do iatom=1, natom
          tempchar=atom_symbols(iatom)
          write(992,'(a8, 2x, 3(f12.6,3x))') tempchar(1:8)
     $      , (tau(ii,iatom)*par*bohr, ii=1,3)
        enddo
        write(992,'(a)') "end atoms_cart"
        
        write(992,'(a, i5)') "num_wann =",num_wann
        write(992,'(a)') "begin unit_cell_cart"
        write(992,'(a)') "Ang"
        write(992,'(3(f14.8,3x))') (rbas(ii,1)*par*bohr, ii=1,3)
        write(992,'(3(f14.8,3x))') (rbas(ii,2)*par*bohr, ii=1,3)
        write(992,'(3(f14.8,3x))') (rbas(ii,3)*par*bohr, ii=1,3)
        write(992,'(a)') "end unit_cell_cart"
        
        write(992,'(a, 3i5)') "mp_grid =", ndiv
        write(992,'(a)') "begin kpoints"
        do ii=1, nqdiv
          write(992, '(3(f14.8, 3x))') kpt_latt(:, ii)
        enddo
        write(992,'(a)') "end kpoints"
        write(992,'(a, i5)') "num_bands =", num_bands
        
        if (ind_low .gt. 1) then
          write(992,'(a, i5,a, i5,a, i5)') "exclude_bands = 1 -",
     $      ind_low, ',  ', ind_high, '-', nbndf
        elseif (ind_low .eq. 1) then
          write(992,'(a, i5, a, i5)') "exclude_bands = 1, ",
     $      ind_high, '-', nbndf                
        else
          write(992,'(a, i5,a,i5)') "exclude_bands = ",
     $      ind_high, '-', nbndf
        endif

        write(992,'(a, i5)') "num_iter = ", num_iter
        write(992,'(a, i5)') "dis_num_iter = ", dis_num_iter
        write(992,'(a, f12.6)') "dis_win_max = ", dis_win_max
        write(992,'(a, f12.6)') "dis_win_min = ", dis_win_min
        write(992,'(a, f12.6)') "dis_froz_max = ", dis_froz_max
        write(992,'(a, f12.6)') "dis_froz_min = ", dis_froz_min

        if (file_exists) then
          write(992,'(a)') "bands_plot = .true."
          write(992,'(a)') "begin kpoint_path"
          OPEN (1, file = 'kpath')
          do ii=1, nlines
            read(1, '(a)') tempstr
            write(992, '(a)') adjustl(trim(tempstr))
          enddo
          close(1)
          write(992,'(a)') "end kpoint_path"                        
!          close(992)
        endif

      close(992)
      endif
      end
      
      

      subroutine cal_num_bands(ispin)
      use comwann_mod

      implicit none
      integer, intent(in) :: ispin
      integer :: ind(1),ik
      double precision, allocatable :: tempvec(:)

      ind_low=10000
      ind_high=-10000          
      do ik=1, npnt
        allocate(tempvec(n_bnd(ik,ispin)))
        if (ubi_num .ge. 4) then
          tempvec=(e_qp(:,ik,ispin)-chem_pot)*rydberg          
        else
          tempvec=(e_bnd(:,ik,ispin)-chem_pot)*rydberg
        endif
        ind=maxloc(tempvec,mask=tempvec .lt. dis_win_min)
        if (ind(1) .lt. ind_low) then
          ind_low=ind(1)
        endif
        ind=minloc(tempvec,mask=tempvec .gt. dis_win_max)
        if (ind(1) .gt. ind_high) then
          ind_high=ind(1)
        endif
        deallocate(tempvec)
      enddo
      if (irel .eq. 2) then
        ind_low=(ind_low/2)*2
        ind_high=((ind_high+1)/2)*2+1
      endif
      num_bands=ind_high-ind_low-1
      if (maswrk) then
        write(iun,*) 'num_bands', num_bands
        write(iun,*) 'ind_low', ind_low
        write(iun,*) 'ind_high', ind_high
      endif

      end      

c$$$  subroutine cal_num_bands(ispin)
c$$$  use comwann_mod
c$$$  
c$$$  implicit none
c$$$  integer, intent(in) :: ispin
c$$$  integer :: ind(1),ik, counter,iq
c$$$  double precision :: meanval, std
c$$$  double precision, allocatable :: tempvec(:)
c$$$  
c$$$  ind_low=10000
c$$$  ind_high=-10000          
c$$$  do ik=1, npnt
c$$$  allocate(tempvec(n_bnd(ik,ispin)))
c$$$  if (ubi_num .ge. 4) then
c$$$  tempvec=(e_qp(:,ik,ispin)-chem_pot)*rydberg          
c$$$  else
c$$$  tempvec=(e_bnd(:,ik,ispin)-chem_pot)*rydberg
c$$$  endif
c$$$  ind=maxloc(tempvec,mask=tempvec .lt. dis_win_min)
c$$$  if (ind(1) .lt. ind_low) then
c$$$  ind_low=ind(1)
c$$$  endif
c$$$  ind=minloc(tempvec,mask=tempvec .gt. dis_win_max)
c$$$  if (ind(1) .gt. ind_high) then
c$$$  ind_high=ind(1)
c$$$  endif
c$$$  deallocate(tempvec)
c$$$  enddo
c$$$  
c$$$  if (npnt .gt. 1) then
c$$$  do counter=1, ind_high-ind_low-1
c$$$  meanval=0.0d0
c$$$  std=0.0d0
c$$$  do ik=1, nqdiv
c$$$  iq=i_kref(ik)
c$$$  if (ubi_num .ge. 4) then          
c$$$  meanval=meanval
c$$$  $         +(e_qp(ind_low+counter,iq,ispin)-chem_pot)
c$$$  $          *rydberg
c$$$  else
c$$$  meanval=meanval
c$$$  $         +(e_bnd(ind_low+counter,iq,ispin)-chem_pot)
c$$$  $          *rydberg
c$$$  endif
c$$$  enddo
c$$$  meanval=meanval/dble(nqdiv)
c$$$  
c$$$  
c$$$  do ik=1, nqdiv
c$$$  iq=i_kref(ik)
c$$$  if (ubi_num .ge. 4) then          
c$$$  std=std
c$$$  $         +((e_qp(ind_low+counter,iq,ispin)-chem_pot)
c$$$  $          *rydberg-meanval)**2
c$$$  else
c$$$  std=std
c$$$  $         +((e_bnd(ind_low+counter,iq,ispin)-chem_pot)
c$$$  $          *rydberg-meanval)**2
c$$$  endif
c$$$  enddo
c$$$  std=dsqrt(std/dble(nqdiv))
c$$$  
c$$$  if (std .gt. 0.5d0) then
c$$$  exit
c$$$  endif
c$$$  enddo
c$$$  ind_low=ind_low+counter-1
c$$$  endif
c$$$  
c$$$  if (irel .eq. 2) then
c$$$  ind_low=(ind_low/2)*2
c$$$  ind_high=((ind_high+1)/2)*2+1
c$$$  endif
c$$$  num_bands=ind_high-ind_low-1
c$$$  if (maswrk) then
c$$$  write(iun,*) 'num_bands', num_bands
c$$$  write(iun,*) 'ind_low', ind_low
c$$$  write(iun,*) 'ind_high', ind_high
c$$$  endif
c$$$  
c$$$  end

      subroutine read_local_axis

      use comwann_mod

      implicit none
      include 'mpif.h'      

      integer :: iatom, ii, jj
      logical :: file_exists

      if (maswrk) then
        inquire(file='local_axis.dat',exist=file_exists)
        if (file_exists) then
          open(unit=55, file='local_axis.dat')
          do iatom=1, natom
            read(55, *) ii, (atom_xaxis(jj,iatom), jj=1, 3),
     $        (atom_zaxis(jj,iatom), jj=1, 3)
          enddo
          close(55)
        else
          do iatom=1, natom
            atom_xaxis(:,iatom)=(/1.0d0, 0.0d0, 0.0d0/)
            atom_zaxis(:,iatom)=(/0.0d0, 0.0d0, 1.0d0/)            
          enddo          
        endif
      endif

      call mpi_bcast(atom_xaxis,3*natom,
     $  mpi_double_precision,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(atom_zaxis,3*natom,
     $  mpi_double_precision,0,MPI_COMM_WORLD,mpierr)
      end
