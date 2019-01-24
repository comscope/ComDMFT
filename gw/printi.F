      subroutine printi
**********************************************************
*  Print input/output data information.                  *
**********************************************************
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      CHARACTER (LEN = 12) REAL_CLOCK (3)
	character*256 :: a
	integer :: i,DATE_TIME(8)
      character*6 :: format_str
      integer :: nchr
      if(.not.maswrk) return
      CALL DATE_AND_TIME (REAL_CLOCK(1),REAL_CLOCK(2),REAL_CLOCK(3),
     &	                DATE_TIME)
      write(iun,1020)real_clock(1)(1:4),real_clock(1)(5:6)
     &	          ,real_clock(1)(7:8)
1020  format(1x,'DATE : ',a4,'/',a2,'/',a2)
      write(iun,1021)real_clock(2)(1:2),real_clock(2)(3:4)
     &	          ,real_clock(2)(5:10)
1021  format(1x,'TIME : ',a2,':',a2,':',a5)
c  Print control parameters
      write(iun,1030)
1030  format(1x,'<Control parameters>    :')
	open(1,file='ini',status='old')
	do i=1,100000
	  read(1,'(a256)',err=1,end=1)a
          nchr = len(trim(a))
          write(format_str,'("(a",i0.3,")")')nchr
	  write(iun,format_str)trim(a)
	enddo
1	continue
	close(1)
      end
