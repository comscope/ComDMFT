



      subroutine cal_rotmat_cmplx2real(norb, rotmat)
      implicit none
      integer, intent(in) :: norb
      complex*16, intent(out) :: rotmat(norb,norb)

      integer :: iorb

      rotmat=0.0d0

      rotmat(norb/2+1,norb/2+1)=1.0d0       
      do iorb=1, norb/2
        rotmat(iorb+norb/2+1,iorb+norb/2+1)
     $    =1.0d0/dsqrt(2.0d0)*(-1)**iorb
        rotmat(iorb+norb/2+1,-iorb+norb/2+1)
     $    =1.0d0/dsqrt(2.0d0)
        rotmat(-iorb+norb/2+1,iorb+norb/2+1)
     $    =1.0d0/dsqrt(2.0d0)*dcmplx(0.0d0, 1.0d0)*(-1)**(iorb+1)
        rotmat(-iorb+norb/2+1,-iorb+norb/2+1)
     $    =1.0d0/dsqrt(2.0d0)*dcmplx(0.0d0, 1.0d0)
      enddo
      end      





      SUBROUTINE  Swap_dble(a, b)
      IMPLICIT  NONE
      double precision, INTENT(INOUT) :: a, b
      double precision                :: Temp
      
      Temp = a
      a    = b
      b    = Temp
      END SUBROUTINE  Swap_dble
      
      SUBROUTINE  Swap_int(a, b)
      IMPLICIT  NONE
      integer, INTENT(INOUT) :: a, b
      integer                :: Temp
      
      Temp = a
      a    = b
      b    = Temp
      END SUBROUTINE  Swap_int
      

      
c     SUBROUTINE  Sort_dble(x, Size, orderarray)
c     IMPLICIT  NONE
c     INTEGER, INTENT(IN) :: Size      
c     double precision, INTENT(INOUT) :: x(size)
c     INTEGER, INTENT(out) :: orderarray(Size)
c     INTEGER :: i,orderarray_temp(Size),Location(1)
c     
c     do i=1, size
c     orderarray_temp(i)=i
c     orderarray(i)=i
c     enddo
c     
c     DO i = 1, Size-1          ! except for the last
c     Location = minloc(x(i:size))+i-1
c     CALL  Swap_dble(x(i), x(Location(1)))
c     CALL  Swap_int(orderarray_temp(i), orderarray_temp(Location(1)))
c     END DO
c     
c     DO i = 1, Size-1          ! except for the last
c     Location = minloc(orderarray_temp(i:size))+i-1
c     CALL  Swap_int(orderarray_temp(i), orderarray_temp(Location(1)))
c     CALL  Swap_int(orderarray(i), orderarray(Location(1)))
c     END DO
c     END



      DOUBLE PRECISION function iFactorial(j)
      IMPLICIT NONE
      INTEGER, intent(in) :: j
      INTEGER :: i
      DOUBLE PRECISION :: x
      if (j<0) print *,
     $  "iFactorial defined only for non-negative numbers!"
      x=1
      iFactorial = x
      if (j.eq.1) return
      DO i=2,j
        x = x*i
      END DO
      iFactorial = x
      return
      end function iFactorial

      DOUBLE PRECISION function dFactorial(x)
      IMPLICIT NONE
      DOUBLE PRECISION, intent(in) :: x
      DOUBLE PRECISION, PARAMETER :: spi2 = 0.8862269254527579
      DOUBLE PRECISION :: y, r
      r=1
      y=x
      DO WHILE(y.gt.1.0)
        r= r * y
        y= y -1.
      ENDDO
      IF (abs(y-0.5).LT.1e-10) r = r*spi2
      dFactorial = r
      return
      END function dFactorial

      DOUBLE PRECISION function mone(i)
      INTEGER, intent(in) :: i
      mone = 1 - 2*MOD(abs(i),2)
      return
      end function mone

      DOUBLE PRECISION function Ddelta(j1, j2, j)
      IMPLICIT NONE
      DOUBLE PRECISION, intent(in) :: j1, j2, j
!     function calls
      DOUBLE PRECISION :: dFactorial
      Ddelta = sqrt(dFactorial(j1+j2-j)*dFactorial(j1-j2+j)
     $  *dFactorial(-j1+j2+j)/dFactorial(j1+j2+j+1))
      return
      END function Ddelta

      DOUBLE PRECISION function f3j(j1, m1, j2, m2, j3, m3)
c$$$  wigner 3j symbol 
      IMPLICIT NONE
      DOUBLE PRECISION, intent(in) :: j1, j2, j3, m1, m2, m3
      INTEGER            :: tmin, tmax, t
      DOUBLE PRECISION             :: sum, v1, v2, dn
!     function calls
      DOUBLE PRECISION             :: dFactorial
      DOUBLE PRECISION             :: iFactorial
      DOUBLE PRECISION             :: Ddelta
      DOUBLE PRECISION             :: mone
      f3j=0
      IF (abs(m1+m2+m3) .GT. 1e-10) return
      IF (abs(j1-j2)-1e-14 .GT. j3 .OR. j3 .GT. j1+j2+1e-14) return
      if (abs(m1) .GT. j1 .OR. abs(m2) .GT. j2 .OR. abs(m3) .GT. j3)
     $  return
      tmin = INT(max(max(0.0,j2-j3-m1),j1-j3+m2)+1e-14)
      tmax = INT(min(min(j1+j2-j3,j1-m1),j2+m2)+1e-14)
      sum=0
      DO t=tmin, tmax
        v1 = dFactorial(j3-j2+m1+t)*dFactorial(j3-j1-m2+t)
        v2 = dFactorial(j1+j2-j3-t)*dFactorial(j1-m1-t)
     $    *dFactorial(j2+m2-t)
        sum = sum + mone(t)/(iFactorial(t)*v1*v2)
      END DO
      dn = dFactorial(j1+m1)*dFactorial(j1-m1)*dFactorial(j2+m2)
     $  *dFactorial(j2-m2)*dFactorial(j3+m3)*dFactorial(j3-m3)
      f3j = mone(INT(j1-j2-m3))*Ddelta(j1,j2,j3)*sqrt(dn)*sum
      return
      END function f3j


c     DOUBLE PRECISION function f3j_int(j1, m1, j2, m2, j3, m3)
c     c$$$  wigner 3j symbol 
c     IMPLICIT NONE
c     integer, intent(in) :: j1, j2, j3, m1, m2, m3
c     INTEGER            :: tmin, tmax, t
c     DOUBLE PRECISION             :: sum, v1, v2, dn,
c     $  j1d,j2d,j3d,m1d,m2d,m3d
c     !     function calls
c     DOUBLE PRECISION             :: dFactorial
c     DOUBLE PRECISION             :: iFactorial
c     DOUBLE PRECISION             :: Ddelta
c     DOUBLE PRECISION             :: mone
c     f3j_int=0
c     j1d=j1*1.0d0
c     j2d=j2*1.0d0
c     j3d=j3*1.0d0
c     m1d=m1*1.0d0
c     m2d=m2*1.0d0
c     m3d=m3*1.0d0
c     IF (abs(m1d+m2d+m3d) .GT. 1e-10) return
c     IF (abs(j1d-j2d)-1e-14.GT.j3d.OR.j3d.GT.j1d+j2d+1e-14) return
c     if (abs(m1d).GT.j1d.OR.abs(m2d).GT.j2d.OR.abs(m3d).GT.j3d)
c     $  return
c     tmin = INT(max(max(0.0,j2d-j3d-m1d),j1d-j3d+m2d)+1e-14)
c     tmax = INT(min(min(j1d+j2d-j3d,j1d-m1d),j2d+m2d)+1e-14)
c     sum=0
c     DO t=tmin, tmax
c     v1 = dFactorial(j3d-j2d+m1d+t)*dFactorial(j3d-j1d-m2d+t)
c     v2 = dFactorial(j1d+j2d-j3d-t)*dFactorial(j1d-m1d-t)
c     $    *dFactorial(j2d+m2d-t)
c     sum = sum + mone(t)/(iFactorial(t)*v1*v2)
c     END DO
c     dn = dFactorial(j1d+m1d)*dFactorial(j1d-m1d)*dFactorial(j2d+m2d)
c     $  *dFactorial(j2d-m2d)*dFactorial(j3d+m3d)*dFactorial(j3d-m3d)
c     f3j_int = mone(INT(j1d-j2d-m3d))*Ddelta(j1d,j2d,j3d)*sqrt(dn)*sum
c     return
c     END function f3j_int      

c     DOUBLE PRECISION function Ggaunt(l1, m1, l2, m2, l3, m3)
c     IMPLICIT NONE
c     INTEGER, intent(in) :: l1, m1, l2, m2, l3, m3
c     DOUBLE PRECISION, PARAMETER   :: pi = 3.14159265358979d0
c     DOUBLE PRECISION :: l1_, l2_, l3_, mm1_, m2_, m3_, zero
c     !     function calls
c     DOUBLE PRECISION             :: f3j
c     DOUBLE PRECISION             :: mone
c     l1_ = l1;   l2_ = l2;   l3_ = l3
c     mm1_ = -m1; m2_ = m2; m3_ = m3
c     zero = 0
c     !     Calculates <Y_{l1m1}|Y_{l2m2}|Y_{l3m3}>
c     if (l1.LT.0 .OR. l2.LT.0 .OR. l3.LT.0) print *,
c     $  "Quantum number l must be non-negative!"
c     Ggaunt = mone(m1)*sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/(4*pi))
c     $  *f3j(l1_,zero,l2_,zero,l3_,zero)*f3j(l1_,mm1_,l2_,m2_,l3_,m3_)
c     return
c     END function Ggaunt

      DOUBLE PRECISION function ClebschG(j,m,j1,m1,j2,m2)
      IMPLICIT NONE
      DOUBLE PRECISION, intent(in) :: j,m,j1,m1,j2,m2
      INTEGER            :: tmin, tmax, t
      DOUBLE PRECISION             :: sum, v1, v2
!     function calls
      DOUBLE PRECISION, external :: iFactorial, dFactorial,mone,Ddelta
      
      ClebschG = 0
      IF (m1+m2 .NE. m) return
      IF ((j .gt. j1+j2) .or. (j .lt. abs(j1-j2))) return      
      tmin = INT(max(max(0.0,j2-j-m1),j1-j+m2)+1e-14)
      tmax = INT(min(min(j1+j2-j,j1-m1),j2+m2)+1e-14)
      sum=0;
      DO t=tmin, tmax
        v1 = sqrt((2*j+1)*dFactorial(j1+m1)*dFactorial(j1-m1)
     $    *dFactorial(j2+m2)*dFactorial(j2-m2)
     $    *dFactorial(j+m)*dFactorial(j-m))
        v2 = iFactorial(t)*dFactorial(j1+j2-j-t)*dFactorial(j1-m1-t)
     $    *dFactorial(j2+m2-t)*dFactorial(j-j2+m1+t)
     $    *dFactorial(j-j1-m2+t)
        sum = sum + mone(t)*v1/v2
      END DO
      ClebschG = sum*Ddelta(j1,j2,j)
      return
      END function ClebschG
      

c     DOUBLE PRECISION function ClebschG_int
c     $  (j_int,m_int,j1_int,m1_int,j2_int,m2_int)
c     ! be careful on the order of j1m1 and j2m2
c     IMPLICIT NONE
c     integer, intent(in) :: j_int,m_int,j1_int,m1_int,j2_int,m2_int
c     DOUBLE PRECISION :: j,m,j1,m1,j2,m2      
c     INTEGER            :: tmin, tmax, t
c     DOUBLE PRECISION             :: sum, v1, v2
c     !     function calls
c     DOUBLE PRECISION             :: iFactorial
c     DOUBLE PRECISION             :: dFactorial
c     DOUBLE PRECISION             :: mone
c     DOUBLE PRECISION             :: Ddelta
c     
c     j=dble(j_int)
c     m=dble(m_int)
c     j1=dble(j1_int)
c     m1=dble(m1_int)
c     j2=dble(j2_int)
c     m2=dble(m2_int)      
c     
c     ClebschG_int = 0
c     IF (m1+m2 .NE. m) return
c     IF ((j .gt. j1+j2) .or. (j .lt. abs(j1-j2))) return
c     tmin = INT(max(max(0.0,j2-j-m1),j1-j+m2)+1e-14)
c     tmax = INT(min(min(j1+j2-j,j1-m1),j2+m2)+1e-14)
c     sum=0;
c     DO t=tmin, tmax
c     v1 = sqrt((2*j+1)*dFactorial(j1+m1)*dFactorial(j1-m1)
c     $    *dFactorial(j2+m2)*dFactorial(j2-m2)
c     $    *dFactorial(j+m)*dFactorial(j-m))
c     v2 = iFactorial(t)*dFactorial(j1+j2-j-t)*dFactorial(j1-m1-t)
c     $    *dFactorial(j2+m2-t)*dFactorial(j-j2+m1+t)
c     $    *dFactorial(j-j1-m2+t)
c     sum = sum + mone(t)*v1/v2
c     END DO
c     ClebschG_int = sum*Ddelta(j1,j2,j)
c     return
c     END function ClebschG_int      






