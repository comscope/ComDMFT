c-----------------------------------------------------------------------
c
      module input_mod
c
c In this module input data is stored that cannot be processed within
c the subroutine readini (the input reader) because the processing
c depends quantities or chunks of memory that are not available at that
c time. Therefore that input data is parked here until it can be
c processed in the appropriate way.
c
        real*8, allocatable :: symmop_rotations(:,:,:)
        real*8, allocatable :: symmop_translations(:,:)
      end
c
c-----------------------------------------------------------------------
