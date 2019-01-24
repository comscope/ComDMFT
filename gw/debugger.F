C> \brief Block execution until a debugger has been attached
C>
C> This subroutine blocks in an infinite loop. By attaching a debugger
C> and setting the variable OBLOCK to .FALSE. the code can be released.
C> Of course this approach is not suitable for large numbers of cores.
C>
      subroutine debugger()
      use manager_mod
      implicit none
#ifdef DEBUGGER
      integer, external :: nnodes
      integer, external :: nodeid
      integer :: istat
      character*(8) :: value
      value = ' '
      call get_environment_variable("GW_ATTACH_DEBUGGER",value,
     +                              STATUS=istat)
      if (istat.eq.2) then
        value(1:1)='Y'
      endif
      dosplot = (nnodes().le.4.and.value(1:1).eq.'Y')
      if (nodeid().eq.0.and.dosplot) then
        write(iun,*)"WARNING: Blocking waiting for attaching debugger"
        flush(iun)
        write(6,*)"WARNING: Blocking waiting for attaching debugger"
        flush(6)
        write(0,*)"WARNING: Blocking waiting for attaching debugger"
        flush(0)
        do while (dosplot) 
          call sleep(1)
        enddo
        write(iun,*)"WARNING: Execution continuing"
        flush(iun)
        write(6,*)"WARNING: Execution continuing"
        flush(6)
        write(0,*)"WARNING: Execution continuing"
        flush(0)
      endif
      call synch(1)
      dosplot=.FALSE.
#endif
      end subroutine debugger
