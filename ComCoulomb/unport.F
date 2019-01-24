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
      use manager_mod
      use parallel_mod
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
c$$$      WRITE(179,*)
c$$$      WRITE(179,'(a)')STRING                           
c$$$      call cpu_time(timeval)
c$$$      WRITE(179,'(a, f20.5)')'time=',timeval
c$$$      WRITE(179,*)
c$$$      call flush(179)
      END                                                         

c      SUBROUTINE log_eachcore(STRING)                                          
c      use manager_mod
c      use parallel_mod
c      implicit none
C     
c      character(len=*),intent(in) :: STRING
c      double precision :: timeval
c
c      WRITE(179,'(a)')trim(STRING)
c      call cpu_time(timeval)
c      WRITE(179,'(a, f20.5)')'time=',timeval
c      WRITE(179,*)
c      call flush(179)
c      END
c

      SUBROUTINE all_log(STRING)                                          
      use manager_mod
      use parallel_mod
      implicit none
C     
      character(len=*),intent(in) :: STRING

c      call log_eachcore(string)
      call timel(string)      
      
      END      

      
c$$$
c$$$
c$$$
c$$$      
c$$$      SUBROUTINE TSECND(TIM)                                            
c$$$C-----------------------------------------------------------------------
c$$$C     ----- THIS ROUTINE PERFORMS CPU AND WALL CLOCK TIMING -----     
c$$$C     
c$$$C     THIS ROUTINE SHOULD SET 'CPU' AND 'WALL' VARIABLES              
c$$$C     TO THE TOTAL ELAPSED CPU AND WALL CLOCK TIMES,                  
c$$$C     MEASURED IN SECONDS.  IN ADDITION, THE CALLING                  
c$$$C     ARGUMENT 'TIM' SHOULD BE SET EQUAL TO 'CPU'.                    
c$$$C     
c$$$C     ON THE FIRST ENTRY, 'CPU0' AND 'WALL0' SHOULD BE SET            
c$$$C     TO THE APPROPRIATE BASE VALUE ON JOB START.                     
c$$$C     
c$$$C     
c$$$C-----------------------------------------------------------------------
c$$$C     * * * *  PC-UNIX SPECIFIC PORTION  * * * *                  
c$$$C     
c$$$      implicit none
c$$$      DOUBLE PRECISION CPU0,WALL0,CPU,WALL,TIM                          
c$$$      integer time                                                      
c$$$      EXTERNAL ETIME                                                    
c$$$      DOUBLE PRECISION ETIME                                            
c$$$      REAL TARRAY(2)                                                    
c$$$      LOGICAL FIRST                                                     
c$$$      COMMON /TIMING/ CPU,WALL                                          
c$$$      SAVE FIRST,CPU0,WALL0                                             
c$$$      DATA FIRST/.TRUE./                                                
c$$$C     
c$$$      IF(FIRST) THEN                                                    
c$$$        FIRST=.FALSE.                                                  
c$$$        WALL0 = TIME()                                                 
c$$$        DUMMY = ETIME(TARRAY)                                          
c$$$        CPU0 = TARRAY(1)+TARRAY(2)                                     
c$$$      END IF                                                            
c$$$C     
c$$$      WALL = TIME()                                                     
c$$$      DUMMY = ETIME(TARRAY)                                             
c$$$      CPU = TARRAY(1)+TARRAY(2)                                         
c$$$      CPU = CPU - CPU0                                                  
c$$$      WALL = WALL - WALL0                                               
c$$$      TIM = CPU                                                         
c$$$      RETURN                                                            
c$$$      end
c$$$
c$$$
c$$$C-----------------------------------------------------------------------
c$$$      SUBROUTINE FLSHBF(LUNIT)                                          
c$$$      use parallel_mod
c$$$C     
c$$$C     
c$$$C     ----- FLUSH THE BUFFER FOR LOGICAL UNIT LUNIT -----            
c$$$C     THIS ROUTINE IS MOST IMPORTANT FOR UNIX SYSTEMS,               
c$$$C     WHERE OUTPUT OTHERWISE STAYS IN BUFFERS FOREVER.               
c$$$C     
c$$$C     NOTE: ON AIX SYSTEMS WITH XLF VERSION 2, THE FLUSH CALL DOES      
c$$$C     NOT WORK WELL, AND SHOULD BE COMMENTED OUT.  THE MOST       
c$$$C     PROMINENT SIGN OF A BAD "FLUSH_" IS GARBAGE CHARACTERS      
c$$$C     IN THE PUNCH FILE.  E.G., CHANGE TO "C---*AIX ..."          
c$$$C     
c$$$      IF (MASWRK) THEN                                                  
c$$$*     AIX  CALL FLUSH_(LUNIT)                                                
c$$$*     APO  CALL FLUSH(LUNIT)                                                 
c$$$*     ARD  CALL FLUSH(LUNIT)                                                 
c$$$*     CEL  CALL FLUSH(LUNIT)                                                 
c$$$*     CRY  CALL FLUSH(LUNIT)                                                 
c$$$*     CVX  CALL FLUSH(LUNIT)                                                 
c$$$*     DEC  CALL FLUSH(LUNIT)                                                 
c$$$        CALL FLUSH(LUNIT)                                                 
c$$$*     INT  CALL FORFLUSH(LUNIT)                                              
c$$$*     KSR  CALL FLUSH(LUNIT)                                                 
c$$$*     SX3  CALL FLUSH(LUNIT)                                                 
c$$$*     SGI  CALL FLUSH(LUNIT)                                                 
c$$$*     STR  CALL FLUSH(LUNIT)                                                 
c$$$*     SUN  CALL FLUSH(LUNIT)                                                 
c$$$*     T3D  CALL FLUSH(LUNIT)                                                 
c$$$*     TMC  CALL FLUSH(LUNIT)                                                 
c$$$        
c$$$      END IF                                                            
c$$$      END

c$$$      subroutine stopgreen()
c$$$      call ending
c$$$      end
