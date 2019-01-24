      SUBROUTINE TIMEL(STRING)                                          
	use manager_mod
	use parallel_mod
      IMPLICIT none                            
      character*36, intent(in) :: STRING
	real*8 :: wall                                   
      CALL TSECND(wall)                                           
      IF (MASWRK) then
	  WRITE(IUN,9008)STRING,' REAL TIME IS ',wall
        FLUSH(iun)
	endif
      RETURN
 9008 FORMAT(a36,a14,F13.3)
      END
      
      
      
      
      
      
      SUBROUTINE TSECND(wall)
      implicit none
      real*8, intent(out) :: wall
	integer :: iff,count,count_rate,count_max
	real*8 :: wall0
      save :: wall0                                            
      DATA iff/0/
      call system_clock(count,count_rate,count_max) 
      wall=real(count)/count_rate
      IF(iff==0) wall0=wall
      iff=1                                                    
      WALL=wall-wall0                           
	end
