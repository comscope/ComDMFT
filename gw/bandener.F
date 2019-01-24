      SUBROUTINE BANDENER(EBND,WEIG,nband)
*************************************************************
*  Calculate one-electron contribution to the total energy  *
*************************************************************
	use manager_mod
	use solid_mod
      IMPLICIT none
	integer :: nband,ispin,k,ibnd
      real*8 :: ebnd,WEIG(nband,NPnT,NSPIN)
	real*8, allocatable :: tos(:,:)
	allocate(tos(0:1,nspin))
      do ispin=1,nspin
        TOS(0,ISPIN)=0.D0
        TOS(1,ISPIN)=0.D0
C   !!! get partial numbers of valence states
        do k=1,npnt                       !!! over k-points
          do ibnd=1,nband                    !!! over bands
            IF(IBND.GE.NBNDE) THEN
              TOS(0,ispin)=TOS(0,ispin)+WEIG(IBND,k,ispin)
              TOS(1,ispin)=TOS(1,ispin)+e_bnd(ibnd,k,ispin)
     &			                     *WEIG(IBND,k,ispin)
            ELSEIF(IBND.LT.NBNDE) THEN
              TOS(0,ispin)=TOS(0,ispin)+WGT(k)/DBLE(NSPIN)
              TOS(1,ispin)=TOS(1,ispin)+e_bnd(ibnd,k,ispin)
     &			                     *WGT(k)/DBLE(NSPIN)
            ENDIF
	    enddo
	  enddo
      enddo
C  PRINT PARTIAL NUMBERS AND DENSITIES OF STATES
      EBND=(TOS(1,1)+TOS(1,NSPIN))/2.D0*DFLOAT(NSPIN)
	deallocate(tos)
      END