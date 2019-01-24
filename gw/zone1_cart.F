	subroutine zone1_cart(q,bas,q1)
	use solid_mod
	implicit none
	real*8, intent(in) :: q(3),bas(3,3)
	real*8, intent(out) :: q1(3)
	real*8 :: det3,dt
	dt=det3(bas(1,1),bas(1,2),bas(1,3))
	q1(1)=det3(q,bas(1,2),bas(1,3))/dt
	q1(1)=mod(q1(1),1.d0)
	if(q1(1)<0.d0) q1(1)=q1(1)+1.d0
	q1(2)=det3(bas(1,1),q,bas(1,3))/dt
	q1(2)=mod(q1(2),1.d0)
	if(q1(2)<0.d0) q1(2)=q1(2)+1.d0
	q1(3)=det3(bas(1,1),bas(1,2),q)/dt
	q1(3)=mod(q1(3),1.d0)
	if(q1(3)<0.d0) q1(3)=q1(3)+1.d0
	q1=matmul(bas,q1)
	end