      subroutine read_ro
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
#ifdef MPI
	include 'mpif.h'
#endif	
      integer :: i_len
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	if(maswrk) then
	  i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_ro',form='unformatted')
        read(3)ro
	  read(3)rointr
        if(magn.eq.2) then
	    read(3)spmt
          read(3)spintr
        endif   !!! only for magn = 2
1       close(3)
	endif
#ifdef MPI
	if(goparr) then
        call brdcst(msgdbl,ro,8*maxmt,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,rointr,16*nplwro*nspin,master,MPI_COMM_WORLD)
        if(magn.eq.2) then
          call brdcst(msgdbl,spmt,8*maxmtb,master,MPI_COMM_WORLD)
          call brdcst(msgdbl,spintr,48*nplwro,master,MPI_COMM_WORLD)
        endif   !!! only for magn = 2
	endif
#endif	
      end