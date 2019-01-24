	subroutine fft1(n1,nx,f,key)
c	key=1  - from K to R
c	key=-1 - from R to K
#ifdef FFTMKL
      use MKL_DFTI
#endif      
	use units_mod
	implicit none
	integer, intent(in) :: n1,nx,key
#ifdef FFTMKL
	complex*16, intent(inout) :: f(n1,nx)
#else	
	real*8, intent(inout) :: f(2*n1,nx)
#endif	
	integer :: j
#ifdef FFTMKL
	integer :: status
      type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle
      complex*16, allocatable :: t(:)
      allocate(t(n1))
c --------------------------------------------------------
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &                              DFTI_COMPLEX,1,n1)
      Status = DftiCommitDescriptor( Desc_Handle )
      if(key==-1) then   !! Forward transform
	  do j=1,nx
	    t=f(:,j)
          Status = DftiComputeForward( Desc_Handle,t)
	    f(:,j)=t
	  enddo
      else if(key==1) then   !! Backward transform
	  do j=1,nx
	    t=f(:,j)
          Status = DftiComputeBackward( Desc_Handle,t)
	    f(:,j)=t
	  enddo
	endif
      Status = DftiFreeDescriptor(Desc_Handle)
      deallocate(t)
#else
	do j=1,nx
	  call cft(f(1,j),f(2,j),n1,n1,n1,2*key)
	enddo
#endif
	end
	
	
	subroutine fft2(n1,n2,nx,f,key)
c	key=1  - from K to R
c	key=-1 - from R to K
#ifdef FFTMKL
      use MKL_DFTI
#endif      
	use units_mod
	implicit none
	integer, intent(in) :: n1,n2,nx,key
#ifdef FFTMKL
	complex*16, intent(inout) :: f(n1*n2,nx)
#else	
	real*8, intent(inout) :: f(2*n1*n2,nx)
#endif	
	integer :: j
#ifdef FFTMKL
	integer :: status,lengths(2)
      type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle
      complex*16, allocatable :: t(:)
      allocate(t(n1*n2))
c --------------------------------------------------------
      lengths(1) = n1
      lengths(2) = n2
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &                              DFTI_COMPLEX,2,lengths)
      Status = DftiCommitDescriptor( Desc_Handle )
      if(key==-1) then   !! Forward transform
	  do j=1,nx
	    t=f(:,j)
          Status = DftiComputeForward( Desc_Handle,t)
	    f(:,j)=t
	  enddo
      else if(key==1) then   !! Backward transform
	  do j=1,nx
	    t=f(:,j)
          Status = DftiComputeBackward( Desc_Handle,t)
	    f(:,j)=t
	  enddo
	endif
      Status = DftiFreeDescriptor(Desc_Handle)
      deallocate(t)
#else
	do j=1,nx
	  call cft(f(1,j),f(2,j),n1*n2,n1,n1,2*key)
	  call cft(f(1,j),f(2,j),n1*n2,n2,n1*n2,2*key)
	enddo
#endif      
	end
	
	
	subroutine fft3(n1,n2,n3,nx,f,key)
c	key=1  - from K to R
c	key=-1 - from R to K
#ifdef FFTMKL
      use MKL_DFTI
#endif      
	use units_mod
	implicit none
	integer, intent(in) :: n1,n2,n3,nx,key
#ifdef FFTMKL
	complex*16, intent(inout) :: f(n1*n2*n3,nx)
#else
	real*8, intent(inout) :: f(2*n1*n2*n3,nx)
#endif	
	integer :: j
#ifdef FFTMKL
	integer :: status,lengths(3)
      type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle
      complex*16, allocatable :: t(:)
      allocate(t(n1*n2*n3))
c --------------------------------------------------------
      lengths(1) = n1
      lengths(2) = n2
      lengths(3) = n3
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &                              DFTI_COMPLEX,3,lengths)
      Status = DftiCommitDescriptor( Desc_Handle )
      if(key==-1) then   !! Forward transform
	  do j=1,nx
	    t=f(:,j)
          Status = DftiComputeForward( Desc_Handle,t)
	    f(:,j)=t
	  enddo
      else if(key==1) then   !! Backward transform
	  do j=1,nx
	    t=f(:,j)
          Status = DftiComputeBackward( Desc_Handle,t)
	    f(:,j)=t
	  enddo
	endif
      Status = DftiFreeDescriptor(Desc_Handle)
      deallocate(t)
#else	
      do j=1,nx
	  call cft(f(1,j),f(2,j),n1*n2*n3,n1,n1,2*key)
	  call cft(f(1,j),f(2,j),n1*n2*n3,n2,n1*n2,2*key)
	  call cft(f(1,j),f(2,j),n1*n2*n3,n3,n1*n2*n3,2*key)
	enddo
#endif      
	end