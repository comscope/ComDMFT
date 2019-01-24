      subroutine fint(f,x,n,t,p)
c     interpolation of function f(x) in point t.
      implicit none
	integer, intent(in) :: n
      real*8, intent(in) :: f(n),x(n),t
	real*8, intent(out) :: p
	integer :: i,j
	real*8 :: s
      p=0.d0
      do i=1,n
        s=f(i)
        do j=1,i-1
          s=s*(t-x(j))/(x(i)-x(j))
        enddo
        do j=i+1,n
          s=s*(t-x(j))/(x(i)-x(j))
        enddo
        p=p+s
      enddo
      end

C>
C> \brief Check if fint will succeed
C>
C> If the values in x are very close together fint might fail with a
C> floating point overflow. This routine returns .false. if this problem
C> is likely to happen, and .true. otherwise.
C>
C> The approach is simple: The largest factor depends on the largest
C> value of (t-x(i)) and the smallest value of (x(i)-x(j)). An upper
C> limit is given by taking the fraction of these two numbers to the
C> (n-1)th power and multiplying it by n. If this estimate exceeds
C> the floating number range then fint cannot be invoked safely.
C>
      logical function ok_fint(f,x,n,t)
      implicit none
      integer, intent(in) :: n    ! the dimension of the arrays
      real*8,  intent(in) :: f(n) ! the function values
      real*8,  intent(in) :: x(n) ! the input coordinates
      real*8,  intent(in) :: t    ! the interpolation point
      real*8 enum  ! the enumerator
      real*8 denom ! the deminator
      real*8 powr  ! the power for the n-th order root
      integer i,j  ! counters
c
      enum  = 0.0d0
      denom = huge(0.0d0)
c
      do i = 1, n
        enum = max(enum,abs(t-x(i)))
        do j = 1, i-1
          denom = min(denom,abs(x(i)-x(j)))
        enddo
        do j = i+1, n
          denom = min(denom,abs(x(i)-x(j)))
        enddo
      enddo
c
      powr = 1.0d0/dble(n-1)
      ok_fint = (n*enum.lt.(denom*(huge(0.0d0)**powr)))
c
      return
      end
