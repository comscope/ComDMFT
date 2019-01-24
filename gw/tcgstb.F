      subroutine main_group(world_group)
#ifdef MPI
      include 'mpif.h'
#endif
      integer world_group
C
C        RETURN THE index of group associated with MPI_WORLD_COMM
C
      world_group = 1
C
#ifdef MPI
      CALL MPI_COMM_GROUP(MPI_COMM_WORLD,world_group,IERROR)
#endif      
      END
      
      
      
      
      
	subroutine group_incl(group_in,ndim,irank,group_out)
#ifdef MPI
      include 'mpif.h'
#endif
      integer, intent(in) :: group_in,ndim,irank(ndim)
      integer, intent(out) :: group_out
C
C        RETURN THE index of group included in group_in
C
      group_out = 1
C
C
#ifdef MPI
      CALL MPI_Group_incl(group_in,ndim,irank,group_out,ierror)
#endif      
      END
	
	
	
	
	subroutine comm_create(group_in,comm_out)
#ifdef MPI
      include 'mpif.h'
#endif
      integer, intent(in) :: group_in
      integer, intent(out) :: comm_out
C
C        Creates the communicator for the given group inside 
!        the MPI_Comm_world
!
      comm_out = 1
C
C
#ifdef MPI
      CALL MPI_COMM_CREATE(MPI_COMM_WORLD,group_in,comm_out,IERROR)
#endif      
      END
	
	
	
	
	
	subroutine cart_create(comm_in,ndims,dims,comm_cart)
#ifdef MPI
      include 'mpif.h'
#endif
      integer, intent(in) :: comm_in,ndims,dims(ndims)
      integer, intent(out) :: comm_cart
	logical periods(ndims),reorder(ndims)
C
C        Creates the cartesian communicator inside 
!        the communicator comm_in
C
      comm_cart = 1
C
C
	periods=.false.
	reorder=.false.
#ifdef MPI
      CALL MPI_Cart_CREATE(comm_in,ndims,dims,periods,reorder,comm_cart,
     &                     IERROR)
#endif     
      END
      
      
      
      
	subroutine cart_sub(comm_cart,remained_dims,comm_cart_new)
#ifdef MPI
      include 'mpif.h'
#endif
      integer, intent(in) :: comm_cart
	logical, intent(in) :: remained_dims(*)
      integer, intent(out) :: comm_cart_new
C
C        Partitions communicator into subgroups with lower dimensions
C
      comm_cart_new = 1
C
C
#ifdef MPI
      CALL MPI_Cart_sub(comm_cart,remained_dims,comm_cart_new,IERROR)
#endif      
      END
      
      
      
      
	subroutine cart_coords(comm_cart,rank_cart,ndims,coord1,coord2)
#ifdef MPI
      include 'mpif.h'
#endif
      integer, intent(in) :: comm_cart,rank_cart,ndims
      integer, intent(out) :: coord1,coord2
	integer :: coords(2)
C
C        Gives the cartesian coordinates of the process with RANK
C
      coord1 = 1
	coord2 = 1
C
C
#ifdef MPI
      CALL MPI_Cart_coords(comm_cart,rank_cart,coords,IERROR)
      coord1 = coords(1)
      coord2 = coords(2)
#endif      
      END
      
      
      
      
	subroutine comm_rank(comm_in,rank)
#ifdef MPI
      include 'mpif.h'
#endif
      integer, intent(in) :: comm_in
      integer, intent(out) :: rank
C
C        Gives the rank of the process in communicator
C
      rank = 0
C
C
#ifdef MPI
      CALL MPI_Comm_rank(comm_in,rank,IERROR)
#endif      
      END
      
      
      
      
	subroutine cart_rank(comm_in,coord,rank)
#ifdef MPI
      include 'mpif.h'
#endif
      integer, intent(in) :: comm_in,coord
      integer, intent(out) :: rank
C
C        Gives the rank of the process in communicator
C
      rank = 0
C
C
#ifdef MPI
      CALL MPI_Cart_rank(comm_in,coord,rank,IERROR)
#endif      
      END
      
      
      
      
	
	subroutine isend(buf,len,idest,comm)
#ifdef MPI
      include 'mpif.h'
#endif
      integer, intent(in) :: len,idest,buf(len),comm
#ifdef MPI
      CALL MPI_SEND(BUF,LEN,MPI_INTEGER,idest,0,comm,IERROR)
#endif      
      END
	
	
	
	
	subroutine irecv(buf,len,ifrom,comm)
#ifdef MPI
      include 'mpif.h'
#endif
      integer, intent(in) :: len,ifrom,buf(len),comm
#ifdef MPI
      CALL MPI_RECV(BUF,LEN,MPI_INTEGER,ifrom,0,comm,istatus,IERROR)
#endif      
      END
	
	
	
	
	
	SUBROUTINE BRDCST(ITYPE,BUF,LENBUF,IFROM,comm)
	integer :: comm,itype,lenbuf,ifrom
      DOUBLE PRECISION BUF(*)
#ifdef MPI
C
C        BROADCAST -LENBUF- BYTES FROM -BUF- ON NODE -IFROM-
C        TO the NODES, belonging to communicator COMM.
C
      INCLUDE 'mpif.h'
      CHARACTER*3 MSGTYP
#ifdef DEBUG
      integer local_len
      local_len = lenbuf
      CALL MPI_BCAST(local_len,1,MPI_INTEGER,0,
     &               comm,IERROR)
      if (local_len.ne.lenbuf) then
        write(*,*)'BRDCST: local and remote dimensions mismatch',
     &            lenbuf,local_len
        flush(6)
        call MPI_Abort(comm,10,ierror)
      endif
#endif
      IF(ITYPE.GE. 65536) MSGTYP='DBL'
      IF(ITYPE.GE.131072) MSGTYP='INT'
      IF(ITYPE.GE.262144) MSGTYP='CHR'
      IF(MSGTYP.EQ.'DBL') THEN
         LEN = LENBUF/8
         CALL MPI_BCAST(BUF,LEN,MPI_DOUBLE_PRECISION,IFROM,
     *                  comm,IERROR)
      END IF
      IF(MSGTYP.EQ.'INT') THEN
         LEN = LENBUF/4
         CALL MPI_BCAST(BUF,LEN,MPI_INTEGER,IFROM,
     *                  comm,IERROR)
      END IF
      IF(MSGTYP.EQ.'CHR') THEN
         CALL MPI_BCAST(BUF,LENBUF,MPI_CHARACTER,IFROM,
     *                  comm,IERROR)
      END IF
#endif      
      END
	
	
	
	
	
	
	
	
	SUBROUTINE BRDCSTc(BUF,LENBUF,IFROM,comm)
	use parallel_mod
#ifdef MPI
      include 'mpif.h'
#endif
	integer :: comm,lenbuf,ifrom
      character*1 BUF(lenbuf)
#ifdef MPI
	integer, allocatable :: ich(:)
	allocate(ich(lenbuf))
C
C        BROADCAST -LENBUF- BYTES FROM -BUF- ON NODE -IFROM-
C        TO ALL OTHERS NODES, SYNCHRONIZING ON MESSAGE RECEIPT.
C
	if(me==ifrom) ich=iachar(buf)
      CALL MPI_BCAST(ich,LENbuf,MPI_INTEGER,IFROM,comm,IERROR)
	buf=achar(ich)
	deallocate(ich)
#endif	
      END
      
      
      
      
      
      
      SUBROUTINE DGOP(X,LENX,OP,comm)
      use parallel_mod
      CHARACTER*3, intent(in) :: OP
      integer, intent(in) :: comm,lenx
      DOUBLE PRECISION, intent(inout) :: X(LENX)
#ifdef MPI
C
C     PERFORM A GLOBAL OPERATION -OP- ON THE ARRAY -X- OF LENGTH -LENX-
C
      INCLUDE 'mpif.h'
#ifdef DEBUG
      integer local_len
      local_len = lenx
      CALL MPI_BCAST(local_len,1,MPI_INTEGER,0,
     &               comm,IERROR)
      if (local_len.ne.lenx) then
        write(*,*)'DGOP: local and remote dimensions mismatch',
     &            lenx,local_len
        flush(6)
        call MPI_Abort(comm,20,ierror)
      endif
#endif
      IF(OP(3:3).EQ.'+') THEN
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,x,lenx,MPI_DOUBLE_PRECISION,
     &                     MPI_SUM,comm,ierror)
      else IF(OP(1:3).EQ.'min') THEN
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,X,1,MPI_DOUBLE_PRECISION,
     &                     MPI_MIN,comm,IERROR)
      else IF(OP(1:3).EQ.'max') THEN
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,X,1,MPI_DOUBLE_PRECISION,
     &                     MPI_MAX,comm,IERROR)
      ELSE
        IF (MASWRK) THEN
          WRITE(6,*) 'DGOP/TCGSTB DOES NOT SUPPORT GLOBAL OP',OP
        END IF
        CALL MPI_ABORT
      END IF
#endif      
      END
      
      
      
      
      
      
      
      
      
      SUBROUTINE IGOP(X,LENX,OP,comm)
      use parallel_mod
      CHARACTER*3, intent(in) :: OP
      integer, intent(in) :: comm,lenx
      integer, intent(inout) :: X(LENX)
#ifdef MPI
      INCLUDE 'mpif.h'
#ifdef DEBUG
      integer local_len
      local_len = lenx
      CALL MPI_BCAST(local_len,1,MPI_INTEGER,0,
     &               comm,IERROR)
      if (local_len.ne.lenx) then
        write(*,*)'IGOP: local and remote dimensions mismatch',
     &            lenx,local_len
        flush(6)
        call MPI_Abort(comm,20,ierror)
      endif
#endif
      IF(OP(3:3).EQ.'+') THEN
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,x,lenx,MPI_INTEGER,MPI_SUM,comm,
     &                     IERROR)
      else IF(OP(1:3).EQ.'min') THEN
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,X,1,MPI_INTEGER,MPI_MIN,comm,
     &                     IERROR)
      else IF(OP(1:3).EQ.'max') THEN
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,X,1,MPI_INTEGER,MPI_MAX,comm,
     &                     IERROR)
      ELSE
        IF (MASWRK) THEN
          WRITE(6,*) 'IGOP/TCGSTB DOES NOT SUPPORT GLOBAL OP',OP
        END IF
        CALL MPI_ABORT
      END IF
#endif      
      END
      
      
      
      
      
      
      
      
      
      
      INTEGER FUNCTION NNODES()
#ifdef MPI
      include 'mpif.h'
#endif
C
C        RETURN THE TOTAL NUMBER OF NODES RUNNING
C
      NNODES = 1
C
C
#ifdef MPI
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERROR)
      NNODES = ISIZE
#endif      
      END
      
      
      
      
      
      
      INTEGER FUNCTION NODEID()
#ifdef MPI
      include 'mpif.h'
#endif
C
C        RETURN THE NODE IDENTIFICATION NUMBER
C        (0 .LE. NODEID .LE. NNODES-1)
C
      NODEID = 0
C
#ifdef MPI
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERROR)
      NODEID = IRANK
#endif      
      END
      
      
      
      
      
      SUBROUTINE PBEGINF
C
C        INITIALIZE THE PARALLEL ENVIRONMENT AT JOB START
C
#ifdef MPI
      CALL MPI_INIT(IERROR)
#endif      
      END
      
      
      
      
      SUBROUTINE PEND
C
C        TERMINATE THE PARALLEL ENVIRONMENT AT JOB END
C
#ifdef MPI
      CALL MPI_FINALIZE(IERROR)
#endif
      END
      
      
      
      
      SUBROUTINE SYNCH(INUM)
      integer :: inum
C
C        SYNCHRONIZE ALL PROCESS THREADS
C
#ifdef MPI
      INCLUDE 'mpif.h'
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
#endif
      END
      
      
      
      
      
      SUBROUTINE MSGSET(MSGDBL,MSGINT,MSGCHR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR
C
C     ----- SELECT USE OF EXTERNAL DATA REPRESENTATION -----
C     TCGMSG CONVERTS MACHINE FORMAT NUMBERS TO A VENDOR INDEPENDENT
C     TYPE KNOWN AS -XDR- WITHIN TCP/IP.  USING THE POWER OF TWO MASKS
C     PERMITS MESSAGES TO EXCHANGE BETWEEN DIFFERENT VENDOR HARDWARE.
C     THE TRANSLATION TO XDR FORMAT IS COSTLY, AND IS AVOIDED BY USING
C     0S WHENEVER ALL NODES IN THE PARALLEL MACHINE ARE THE SAME.
C
      IF(IXDR.EQ.1) THEN
         MSGDBL =  65536
         MSGINT = 131072
         MSGCHR = 262144
      ELSE
         MSGDBL = 0
         MSGINT = 0
         MSGCHR = 0
      END IF
#ifdef MPI
      MSGDBL =  65536
      MSGINT = 131072
      MSGCHR = 262144
#endif
      END
      
      
      
      SUBROUTINE ENDING                                                 
C                                                                       
C     ----- TERMINATE EXECUTION SMOOTHLY -----                          
C                                                                       
C     CLEAN UP PARALLEL EXECUTION                                       
C                                                                       
      CALL SYNCH(10)                                                    
      CALL PEND()                                                       
      STOP
C                                                           
      END
      
      
      
      SUBROUTINE flshbf(lunit)
      implicit none
      integer, intent(in) :: lunit
      call flush(lunit)
      END            
