c-----------------------------------------------------------------------
!>
!> \brief Module to read input data from stdin in free format.
!>
!> The issue is to be able to read free format input while maintaining
!> compatibility with the original fixed format input. In particular
!> inputs like:
!> \verbatim
!>   admix=0.050
!> \endverbatim
!> cause problems. In fixed format the positions of the floating point
!> number allow reading the data correctly, but in free format the whole
!> line is one token (and not in the correct format for floating point
!> numbers). Hence we need the ability to read a line and tokenize it
!> using separators "=", ",", "(", ")" and " ". Subsequently we need to
!> be able to read individual tokens of appropriate types.
!>
      module inp_mod
!>
!>      Maximum line length
        integer, parameter, private :: maxchr = 256
!>
!>      Maximum number of tokens on a line
        integer, parameter, private :: maxtok = 50
!>
!>      The line format string
        character*6, parameter, private :: fmtline = '(a256)'
!>
!>      The current line
        character*(maxchr), private :: line
!>
!>      The number of tokens on the line
        integer, private :: numtok
!>
!>      The number of the current token
        integer, private :: curtok
!>
!>      Token start indeces
        integer, private :: itokbeg(maxtok)
!>
!>      Token end indeces
        integer, private :: itokend(maxtok)
!>
!>      If .true. all errors a fatal, otherwise return an error
        logical, private :: ohardfail = .false.
!>
!>      The unit number for error messages
        integer, private :: iuerr = 6
c
c       Error codes
c
!>      Everything is fine
        integer, parameter, public :: ERR_OK   = 0
!>
!>      End-of-file reached
        integer, parameter, public :: ERR_EOF  = 1
!>
!>      End-of-line reached
        integer, parameter, public :: ERR_EOL  = 2
!>
!>      The token is not of the expected data type
        integer, parameter, public :: ERR_TYPE = 3
!>
!>      Unknown error
        integer, parameter, public :: ERR_UNK  = 4
c
      contains
!>
!>      \brief Read and tokenize a line of the input file
!>
        subroutine inp_line(iunit,ierr,separators)
        implicit none
!>
!>      The input file unit number
        integer, intent(in)  :: iunit
!>
!>      The error code
        integer, intent(out) :: ierr
!>
!>      A string of token separators
        character*(*), optional, intent(in) :: separators
!
        integer :: length  ! the length of the line in characters
        logical :: intoken ! are we currently in a token?
        integer :: i       ! counter for looping over characters
        character*80 :: seps ! separators
        if (present(separators)) then
          seps = separators
        else
          seps = '= (),'
        endif
        numtok  = 0
        curtok  = 1
        ierr    = ERR_UNK
        intoken = .false.
        if (ohardfail) then
          read(iunit,fmtline,end=100)line
        else
          read(iunit,fmtline,end=100,err=200)line
        endif
        length = len(trim(line))
        do i = 1, length
          if (intoken) then
            if (index(trim(seps),line(i:i)).ne.0) then
              intoken = .false.
              itokend(numtok) = i-1
            endif
          else
            if (index(trim(seps),line(i:i)).eq.0) then
              intoken = .true.
              numtok = numtok + 1
              itokbeg(numtok) = i
            endif
          endif
        enddo
        if (intoken) then
          itokend(numtok) = length
        endif
c       Some numbers to represent having reached the end of the line
c       This helps subroutine inp_print
        itokbeg(numtok+1) = itokend(numtok)+5
        itokend(numtok+1) = itokbeg(numtok+1)+5
        ierr = ERR_OK
        return
 100    ierr = ERR_EOF
 200    return
        end subroutine inp_line
!>
!>      \brief Return the number of tokens on the current line
!>
        integer function inp_len()
        implicit none
        inp_len = numtok
        return
        end function inp_len
!>
!>      \brief Return the current position on the line
!>
!>      Return the current position, in terms of tokens, on the line.
!>
        integer function inp_pos()
        implicit none
        inp_pos = curtok
        return
        end function inp_pos
!>
!>      \brief Read an integer
!>
!>      Read an integer value from the current input line and
!>      move to the next token.
!>
        subroutine inp_int(ival,ierr)
        implicit none
!>
!>      The integer value read
        integer, intent(out) :: ival
!>
!>      The error code
        integer, intent(out) :: ierr
        ierr = ERR_EOL
        if (curtok.gt.numtok) then
          if (ohardfail) then
            call inp_print(iuerr)
            write(iuerr,*)'ERROR: reading beyond the end of the line'
            stop
          endif
          return
        endif
        ierr = ERR_TYPE
        if (ohardfail) then
          read(line(itokbeg(curtok):itokend(curtok)),*,err=200)ival
        else
          read(line(itokbeg(curtok):itokend(curtok)),*,err=100)ival
        endif
        curtok = curtok + 1
        ierr = ERR_OK
 100    return
 200    call inp_print(iuerr)
        write(iuerr,*)'ERROR: Could not read integer'
        stop
        end subroutine inp_int
!>
!>      \brief Read a logical
!>
!>      Read a logical value from the current input line and
!>      move to the next token.
!>
        subroutine inp_log(oval,ierr)
        implicit none
!>
!>      The logical value read
        logical, intent(out) :: oval
!>
!>      The error code
        integer, intent(out) :: ierr
        ierr = ERR_EOL
        if (curtok.gt.numtok) then
          if (ohardfail) then
            call inp_print(iuerr)
            write(iuerr,*)'ERROR: reading beyond the end of the line'
            stop
          endif
          return
        endif
        ierr = ERR_TYPE
        if (ohardfail) then
          read(line(itokbeg(curtok):itokend(curtok)),*,err=200)oval
        else
          read(line(itokbeg(curtok):itokend(curtok)),*,err=100)oval
        endif
        curtok = curtok + 1
        ierr = ERR_OK
 100    return
 200    call inp_print(iuerr)
        write(iuerr,*)'ERROR: Could not read logical'
        stop
        end subroutine inp_log
!>
!>      \brief Read a double precision value
!>
!>      Read an double precision value from the current input line and
!>      move to the next token.
!>
        subroutine inp_dbl(val,ierr)
        implicit none
!>
!>      The double precision value read
        double precision, intent(out) :: val
!>
!>      The error code
        integer, intent(out) :: ierr
        ierr = ERR_EOL
        if (curtok.gt.numtok) then
          if (ohardfail) then
            call inp_print(iuerr)
            write(iuerr,*)'ERROR: reading beyond the end of the line'
            stop
          endif
          return
        endif
        ierr = ERR_TYPE
        if (ohardfail) then
          read(line(itokbeg(curtok):itokend(curtok)),*,err=200)val
        else
          read(line(itokbeg(curtok):itokend(curtok)),*,err=100)val
        endif
        curtok = curtok + 1
        ierr = ERR_OK
 100    return
 200    call inp_print(iuerr)
        write(iuerr,*)'ERROR: Could not read double precision'
        stop
        end subroutine inp_dbl
!>
!>      \brief Read a string value
!>
!>      Read the current token as a string and move to the next token.
!>
        subroutine inp_str(sval,ierr)
        implicit none
!>
!>      The string value read
        character*(*), intent(out) :: sval
!>
!>      The error code
        integer, intent(out) :: ierr
        integer :: i, j
        ierr = ERR_EOL
        if (curtok.gt.numtok) then
          if (ohardfail) then
            call inp_print(iuerr)
            write(iuerr,*)'ERROR: reading beyond the end of the line'
            stop
          endif
          return
        endif
        ierr = ERR_TYPE
        if (itokend(curtok)-itokbeg(curtok)+1.gt.len(sval)) then
          if (ohardfail) then
            go to 200
          else
            go to 100
          endif
        endif
        j = 0
        sval = " "
        do i = itokbeg(curtok), itokend(curtok)
          j = j + 1
          sval(j:j) = line(i:i)
        enddo
        curtok = curtok + 1
        ierr = ERR_OK
 100    return
 200    call inp_print(iuerr)
        write(iuerr,*)'ERROR: Could not read string'
        stop
        end subroutine inp_str
!>
!>      \brief Skip a number of tokens
!>
!>      Skip a number of tokens on the input line. A negative 
!>      value moves the current token back. The lowest token
!>      is 1 and the highest is numtok+1.
!>
        subroutine inp_skip(inum,ierr)
        implicit none
!>
!>      The number of tokens to skip
        integer, intent(in)  :: inum
!>
!>      The error code
        integer, intent(out) :: ierr
        curtok = min(numtok+1,max(1,curtok+inum))
        ierr = ERR_EOL
        if (curtok.gt.numtok) return
        ierr = ERR_OK
        end subroutine inp_skip
!>
!>      \brief Return whether the end of this line was reached
!>
        logical function inp_end()
        implicit none
        inp_end = curtok.gt.numtok
        end function inp_end
!>
!>      \brief Manage the hardfail setting
!>
!>      If oval is .True. turn on hard fail error handling,
!>      otherwise handle input errors gracefully.
!>
        subroutine inp_hardfail(oval)
        implicit none
!>
!>      The requested hardfail setting
        logical, intent(in) :: oval
        ohardfail = oval
        end subroutine inp_hardfail
!>
!>      \brief Manage the error message output unit
!>
        subroutine inp_errunit(iunit)
        implicit none
!>
!>      The requested unit number for error messages
        integer, intent(in) :: iunit
        iuerr = iunit
        end subroutine inp_errunit
!>
!>      \brief Print the current line
!>
!>      Print the current line of input and the current position within
!>      that line.
!>
        subroutine inp_print(iunit)
        implicit none
!>
!>      The file unit number for the output
        integer, intent(in) :: iunit
c
        integer :: ii
c
        write(iunit,'(1x,a)')trim(line)
        write(iunit,'(" ",$)')
        do ii = 1,itokbeg(curtok)-1
          write(iunit,'(a,$)')"-"
        enddo
        do ii = itokbeg(curtok),itokend(curtok)
          write(iunit,'(a,$)')"^"
        enddo
        do ii = itokend(curtok)+1,itokend(numtok)
          write(iunit,'(a,$)')"-"
        enddo
        write(iunit,*)
        end subroutine inp_print

      end module inp_mod
c
c-----------------------------------------------------------------------
