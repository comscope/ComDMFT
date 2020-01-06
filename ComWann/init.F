      subroutine init
*************************************************************
*     control module <init> for reading all the input data   *
*************************************************************
      use comwann_mod
      implicit none
       call read_comwann_ini
#ifdef USE_HDF5
       call read_wan_bni_hdf5
       call read_wan_bni2_hdf5
#else
       call read_wan_bni
       call read_wan_bni2       
#endif
!      call timel('******** ComWann started ***************')
!      call printi
      end
