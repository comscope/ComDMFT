      subroutine init
*************************************************************
*    control module <init> for reading all the input data   *
*************************************************************
	use manager_mod
	use parallel_mod
      implicit none
	call msgset(msgdbl,msgint,msgchr)
      call readini
      call timel('******** QFT started ***************')
      call printi
      end