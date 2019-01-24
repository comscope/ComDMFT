      subroutine beging
	use manager_mod
	use parallel_mod
	implicit none
	integer :: ixdr,nodeid,nnodes        
C                                                                       
C INITIALIZE PARALLEL                                                   
C                                                                       
      CALL PBEGINF()                                                    
      IXDR = 1                                                          
      ME = NODEID()                                                     
      MASTER = 0                                                        
      NPROC = NNODES()                                                  
      IF(NPROC.GT.1) THEN                                               
         GOPARR = .TRUE.                                                
      ELSE                                                              
         GOPARR = .FALSE.                                               
      END IF                                                            
      MASWRK = ME .EQ. MASTER
	end