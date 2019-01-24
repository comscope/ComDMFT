      subroutine goto_end(ifile)
	integer :: i,ifile
      do i=1,1000000
        read(ifile,*,end=1)
      enddo
1     backspace ifile
      end