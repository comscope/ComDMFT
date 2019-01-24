      subroutine integral_pbint_pw(int_pw,v,iq,s)
	use atom_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: iq
      real*8, intent(in) :: v(3)
      complex*16, intent(in) :: s(nplwgw(iq),nplwgw(iq))
      complex*16, intent(out) :: int_pw(nplwgw(iq))
	integer :: iatom,isort,i,i0
      real*8 :: pi4,qq,pi2,phase,v1(3),pi2a,qq3,x,f,v0(3)
      complex*16 :: cc
      complex*16, allocatable :: tmp(:)
	pi2=pi+pi
	pi4=4.d0*pi
	pi2a=pi2/par
	v0=v-pnt(:,iq)
	allocate(tmp(nplwgw(iq)))
	tmp=(0.d0,0.d0)
	do i=1,nplwgw(iq)
	  i0=indpw_gw(i,iq)
	  v1=v0-gbs(:,i0)
	  qq=sqrt(dot_product(v1,v1))*pi2a
	  if(qq<1.d-9) then
	    tmp(i)=amegaint
	  else
	    qq3=pi4/qq**3
	    do isort=1,nsort
	      x=qq*smt(isort)
	      f=qq3*(sin(x)-x*cos(x))
	      do iatom=1,natom
	        if(is(iatom)==isort) then
	          phase=pi2*dot_product(v1,tau(:,iatom))
	          cc=dcmplx(cos(phase),sin(phase))
	          tmp(i)=tmp(i)-cc*f
	        endif
	      enddo
	    enddo
	  endif
	enddo
	call zgemv('n',nplwgw(iq),nplwgw(iq),(1.d0,0.d0),s,nplwgw(iq),tmp,
     &	         1,(0.d0,0.d0),int_pw,1)
	deallocate(tmp)
      end