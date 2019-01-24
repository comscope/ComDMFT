      subroutine radmsh
*************************************************************
*   module prepares radial mesh in each  mt-sphere          *
*************************************************************
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: isort,irad,l
	real*8 :: betta
      do isort=1,nsort
        betta=smt(isort)/(exp(h(isort)*nrad(isort))-1.d0)
*  define radial mesh arrays
        do irad=0,nrad(isort)
          r(irad,isort)=betta*(exp(h(isort)*irad)-1.d0)
          dr(irad,isort)=betta*exp(h(isort)*irad)
          weight(irad,isort)=4.d0*pi*r(irad,isort)**2
     &                                 *dr(irad,isort)
          rl(irad,0,isort)=1.d0
          do l=1,max_rl
            rl(irad,l,isort)=r(irad,isort)**l
          enddo
        enddo
        if(maswrk) then
          write(iun,'(a30,2f11.8)')' Step 1 ; Step N on rad. mesh ',
     &                               r(1,isort),
     &                               r(nrad(isort),isort)
     &                              -r(nrad(isort)-1,isort)
        endif
      enddo
      end