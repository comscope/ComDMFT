      Program QFT
      use comwann_mod
!      use gitinfo_mod

************************************************************
*     The main module                             *
************************************************************

      implicit none
      include 'mpif.h'
      character*300 :: ch
      character*300 :: chbk
!      logical :: file_exist

      call beging
!      mem_using=0.d0
!      mem_max=0.d0
!      mem_where='qft'


c      if (me .lt. 10) then
c        write(ch, '(a, i1, a)') "check_", me, ".dat"
c      elseif (me .lt. 100) then
c        write(ch, '(a, i2, a)') "check_", me, ".dat"
c      elseif (me .lt. 1000) then
c        write(ch, '(a, i3, a)') "check_", me, ".dat"
c      elseif (me .lt. 10000) then
c        write(ch, '(a, i4, a)') "check_", me, ".dat"
c      endif


c      if (me .lt. 10) then
c        write(chbk, '(a, i1, a)') "check_", me, ".bkt"
c      elseif (me .lt. 100) then
c        write(chbk, '(a, i2, a)') "check_", me, ".bkt"
c      elseif (me .lt. 1000) then
c        write(chbk, '(a, i3, a)') "check_", me, ".bkt"
c      elseif (me .lt. 10000) then
c        write(chbk, '(a, i4, a)') "check_", me, ".bkt"
c      endif


c      open(179, file=trim(ch))
c      open(180, file=trim(chbk))
c      write(179,*)'bkk-test',me
      iun=993
      if(maswrk) then
        open(iun,file='OUT')
!        write(iun,'(a)')  "git info"
!        write(iun,'(a)')  gitdate
!        write(iun,'(a)')  gitbranch
        write(iun,*)  
        write(iun,*)  
      endif

      call init
      if(maswrk) write(iun,*)' Number of processors is ',nproc

!      if (me .lt. 10) then
!        write(me_ch, '(i1)') me
!      elseif (me .lt. 100) then
!        write(me_ch, '(i2)') me
!      elseif (me .lt. 1000) then
!        write(me_ch, '(i3)') me
!      elseif (me .lt. 10000) then
!        write(me_ch, '(i4)') me
!      else 
!        write(me_ch, '(i5)') me
!      endif

      call start_par_solid
      call set_memory

#ifdef USE_HDF5
      call read_chemical_potential_rst_hdf5
      call read_wfn_effhf_rst_hdf5
#else
      call read_chemical_potential_rst
      call read_wfn_effhf_rst
#endif
      call wannier_interface
c      close(179)
c      close(180)
      if(maswrk) then
        close(iun)

      endif
      call ending_qft
      end
