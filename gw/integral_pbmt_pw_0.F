      subroutine integral_pbmt_pw_0(mt_pw_0,qq)
	use atom_mod
	use solid_mod
	use units_mod
      implicit none
      real*8, intent(in) :: qq
      complex*16, intent(out) :: mt_pw_0(maxlfpb,nsort)
	integer :: isort,l,ie,i
      real*8 :: pi4,res
      complex*16 :: cc
	pi4=4.d0*pi
	do isort=1,nsort
        do l=0,lmpb(isort)
          cc=(0.d0,1.d0)**l*pi4
	    do ie=1,ntle_pb(l,isort)
	      i=ind_prod0(ie,l,isort)
            call integral_m_j(ie,l,isort,qq,res)
            mt_pw_0(i,isort)=cc*res
	    enddo  !! over ie
	  enddo !! over l
	enddo   !! over isort
      end