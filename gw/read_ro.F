c$$$      subroutine read_ro
c$$$	use atom_mod
c$$$	use manager_mod
c$$$	use parallel_mod
c$$$	use solid_mod
c$$$      implicit none
c$$$#ifdef MPI
c$$$	include 'mpif.h'
c$$$#endif	
c$$$      integer :: i_len
c$$$c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c$$$	if(maswrk) then
c$$$	  i_len=len_trim(allfile)
c$$$        open(3,file=allfile(1:i_len)//'_ro',form='unformatted')
c$$$        read(3)ro
c$$$	  read(3)rointr
c$$$        if(magn.eq.2) then
c$$$	    read(3)spmt
c$$$          read(3)spintr
c$$$        endif   !!! only for magn = 2
c$$$1       close(3)
c$$$	endif
c$$$#ifdef MPI
c$$$	if(goparr) then
c$$$        call brdcst(msgdbl,ro,8*maxmt,master,MPI_COMM_WORLD)
c$$$        call brdcst(msgdbl,rointr,16*nplwro*nspin,master,MPI_COMM_WORLD)
c$$$        if(magn.eq.2) then
c$$$          call brdcst(msgdbl,spmt,8*maxmtb,master,MPI_COMM_WORLD)
c$$$          call brdcst(msgdbl,spintr,48*nplwro,master,MPI_COMM_WORLD)
c$$$        endif   !!! only for magn = 2
c$$$	endif
c$$$#endif	
c$$$      end
