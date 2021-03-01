      subroutine init
*************************************************************
*    control module <init> for reading all the input data   *
*************************************************************
	use manager_mod
	use parallel_mod
      implicit none
      call msgset(msgdbl,msgint,msgchr)
c$$$  create_directory
#ifndef USE_HDF5      
      checkpoint_dir='./checkpoint/'
      if (maswrk) then
        call system('mkdir -p '//adjustl(trim(checkpoint_dir)))
      endif
#endif      
      call readini
      call timel('******** QFT started ***************')
      call printi
      end
