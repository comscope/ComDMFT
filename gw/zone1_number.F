	subroutine zone1_number(q,bas,nd,ind)
c     Usage:
c     For PNT  ---->  RB0,NDIV
c     For RBAS ---->  GBAS,NDIV
c     For r    ---->  Q0B0,NRDIV
c     For GBAS ---->  RBAS,NRDIV
	use solid_mod
	implicit none
	integer, intent(in) :: nd(3)
	real*8, intent(in) :: q(3),bas(3,3)
	integer, intent(out) :: ind
	integer :: i1,i2,i3,ierr
	real*8 :: rv
	ierr=0
      rv=dot_product(q,bas(:,1))
      i1=nint(rv)
      if(abs(rv-dble(i1))>1.d-6) ierr=1
	i1=mod(i1,nd(1))
	if(i1.lt.0) i1=i1+nd(1)
	ind=i1+1
      rv=dot_product(q,bas(:,2))
      i2=nint(rv)
      if(abs(rv-dble(i2))>1.d-6) ierr=1
	i2=mod(i2,nd(2))
	if(i2.lt.0) i2=i2+nd(2)
	ind=ind+nd(1)*i2
      rv=dot_product(q,bas(:,3))
      i3=nint(rv)
      if(abs(rv-dble(i3))>1.d-6) ierr=1
	i3=mod(i3,nd(3))
	if(i3.lt.0) i3=i3+nd(3)
	ind=nd(1)*nd(2)*i3+ind
	if(ierr==1) ind=-ind
	end
	
	subroutine zone1_number_v(q,bas,nd,iv)
	use solid_mod
	implicit none
	integer, intent(in) :: nd(3)
	real*8, intent(in) :: q(3),bas(3,3)
	integer, intent(out) :: iv(3)
	integer :: i1,i2,i3
	real*8 :: ddot
      i1=nint(ddot(3,q,1,bas(1,1),1))
	i1=mod(i1,nd(1))
	if(i1.lt.0) i1=i1+nd(1)
      i2=nint(ddot(3,q,1,bas(1,2),1))
	i2=mod(i2,nd(2))
	if(i2.lt.0) i2=i2+nd(2)
      i3=nint(ddot(3,q,1,bas(1,3),1))
	i3=mod(i3,nd(3))
	if(i3.lt.0) i3=i3+nd(3)
	iv=(/i1,i2,i3/)
	end
	
	subroutine zone1_number2d(q,bas,nd,ind)
c     Usage:
c     For PNT  ---->  RB0,NDIV
c     For RBAS ---->  GBAS,NDIV
c     For r    ---->  Q0B0,NRDIV
c     For GBAS ---->  RBAS,NRDIV
	use solid_mod
	implicit none
	integer, intent(in) :: nd(2)
	real*8, intent(in) :: q(2),bas(2,2)
	integer, intent(out) :: ind
	integer :: i1,i2,ierr
	real*8 :: rv
	ierr=0
      rv=dot_product(q,bas(:,1))
      i1=nint(rv)
      if(abs(rv-dble(i1))>1.d-6) ierr=1
	i1=mod(i1,nd(1))
	if(i1.lt.0) i1=i1+nd(1)
	ind=i1+1
      rv=dot_product(q,bas(:,2))
      i2=nint(rv)
      if(abs(rv-dble(i2))>1.d-6) ierr=1
	i2=mod(i2,nd(2))
	if(i2.lt.0) i2=i2+nd(2)
	ind=ind+nd(1)*i2
	if(ierr==1) ind=-ind
	end