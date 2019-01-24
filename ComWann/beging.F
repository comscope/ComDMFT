      subroutine beging
      use comwann_mod
      implicit none
      include "mpif.h"
      integer :: ierror,irank,isize
C     
C     INITIALIZE PARALLEL                                                   
C     
      CALL MPI_INIT(IERROR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERROR)
      ME = irank
      MASTER = 0                                                        
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERROR)
      NPROC = isize
      IF(NPROC.GT.1) THEN                                               
        GOPARR = .TRUE.                                                
      ELSE                                                              
        GOPARR = .FALSE.                                               
      END IF                                                            
      MASWRK = ME .EQ. MASTER
      end
