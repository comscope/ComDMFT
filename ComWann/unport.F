C-----------------------------------------------------------------------
C     THIS IS THE NASTIEST, TOTALLY UNPORTABLE CODE                     
C     THE VERSIONS WHICH ARE SUPPORTED HERE ARE                         
C     *AIX - ANY IBM RUNNING AIX (SEE IRON.DOC)                    
C     *ALL - ALLIANT FX SERIES (UNIX)                              
C     *AMD - AMDAHL                                                
C     *APO - APOLLO DN 10000 (SEE ALSO *A68 AT BOTTOM)             
C     *ARD - THE ARDENT TITAN PRODUCT LINE FROM STARDENT           
C     *CEL - CELERITY, AND FPS MODEL 500 (UNIX)                    
C     *CRY - CRAY RESEARCH (UNICOS) - SEE ALSO *T3D                
C     *CVX - CONVEX (UNIX)                                         
C     *DEC - DECSTATION/DECSYSTEM (UNIX)                           
C     *DGA - DATA GENERAL AVIION (GREEN HILLS GHF77, ABSOFT F77)   
C     *F2C - F2C/GCC COMPILER ON PC RUNNING UNIX (LINUX,BSD386,ETC)
C     *F77 - A GENERIC, AND MOSTLY DO NOTHING VERSION              
C     *FUJ - FUJITSU UXP/M                                         
C     *HP  - HEWLETT-PACKARD HP/9000 7X0 SERIES                    
C     *IBM - IBM MAINFRAME (VM OR MVS, VS FORTRAN)                 
C     *INT - INTEL PARAGON OR IPSC/860 (UNIX)                      
C     *SX3 - NEC SX-3 (SUPER-UX)                                   
C     *SGI - SILICON GRAPHICS INC.                                 
C     *STR - THE STELLAR PRODUCT LINE FROM STARDENT                
C     *SUN - SUN WORKSTATIONS                                      
C     *T3D - CRAY T3D MACHINES                                     
C     *TMC - THINKING MACHINES (THE CM-5)                          
C     *VAX - VAX/VMS SYSTEMS                                       
C     *DS1 - FTN/LAHEY SYSTEMS                                     
C-----------------------------------------------------------------------
C     
      SUBROUTINE ENDING
      implicit none
      include "mpif.h"
      integer :: ierror

c      write(179,*) 'call ending'
c      call flush(179)

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      CALL MPI_FINALIZE(IERROR)
      stop      
      END                                                               

      SUBROUTINE ENDING_qft
      implicit none
      include "mpif.h"
      integer :: ierror

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      CALL MPI_FINALIZE(IERROR)
      stop            
      END                                                               

      SUBROUTINE TIMEL(STRING)                                          
      use comwann_mod
      implicit none
C     
      character(len=*),intent(in) :: STRING
      double precision :: timeval
      IF (MASWRK) then
        WRITE(IUN,*)
        WRITE(IUN,'(a)')STRING                           
        call cpu_time(timeval)
        WRITE(IUN,'(a, f20.5)')'time=',timeval
        WRITE(IUN,*)
        call flush(iun)
      endif
      END                                                         

      SUBROUTINE log_eachcore(STRING)                                          
      use comwann_mod
      implicit none
C     
      character(len=*),intent(in) :: STRING
      double precision :: timeval

      WRITE(179,'(a)')trim(STRING)
      call cpu_time(timeval)
      WRITE(179,'(a, f20.5)')'time=',timeval
      WRITE(179,*)
      call flush(179)
      END


c      SUBROUTINE all_log(STRING)                                          
c      use comwann_mod
c      implicit none
cC     
c      character(len=*),intent(in) :: STRING
c
c      call log_eachcore(string)
c      call timel(string)      
c      
c      END      
