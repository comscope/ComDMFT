  	subroutine inhomo(ph_in,pi_in,qh_in,qi_in,ph,qh,ppi,qqi,node,
     &	                q,nr,kc,dr,h)
      use manager_mod
      use units_mod
      implicit none
      integer, intent(in) :: nr,kc
      integer, intent(out) :: node
      real*8, intent(in) :: ph_in,pi_in,qh_in,qi_in,ppi(0:nr),qqi(0:nr),
     &                      dr(0:nr),h
      real*8, intent(inout) :: ph(0:nr),qh(0:nr)
      real*8, intent(out) :: q
	integer :: ir
      real*8 :: ph_out,qh_out,pi_out,qi_out,dt,pp,qq,aa,bb,c2,dqdall
      real*8, allocatable :: work(:)
	ph_out=ph(kc)
	qh_out=qh(kc)
	pi_out=ppi(kc)
	qi_out=qqi(kc)
	dt=qh_out*ph_in-ph_out*qh_in
	pp=pi_in-pi_out
	qq=qi_in-qi_out
	aa=(ph_in*qq-pp*qh_in)/dt
	bb=(ph_out*qq-pp*qh_out)/dt
	do ir=0,kc
	  ph(ir)=ppi(ir)+aa*ph(ir)
	  qh(ir)=qqi(ir)+aa*qh(ir)
	enddo
	do ir=kc+1,nr
	  ph(ir)=ppi(ir)+bb*ph(ir)
	  qh(ir)=qqi(ir)+bb*qh(ir)
	enddo
	allocate(work(0:nr))
	do ir=0,nr
	  work(ir)=ph(ir)**2*dr(ir)
      enddo
	if(irel.ge.1) then
	  c2=clight**2
	  do ir=0,nr
	    work(ir)=work(ir)+qh(ir)**2*dr(ir)/c2
        enddo
	endif
      q=dqdall(h,work(0),nr)
	deallocate(work)
	node=0
	do ir=1,nr
	  if(ph(ir-1)*ph(ir).lt.0.d0) node=node+1
      enddo
	end