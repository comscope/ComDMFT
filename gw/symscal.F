      subroutine symscal(fun)
c     nplwro  - number of plane wave for density (in sphere cutoffro)
      use manager_mod
	use solid_mod
	use units_mod
      implicit none
      complex*16, intent(inout) :: fun(nplwro)
	integer :: istar,mm,iv,ig,jv,jj,ii
	real*8 :: x(3),delt,fas,cosn,sinn
	real*8, allocatable :: tmp(:),tmpi(:)
      allocate(tmp(maxstar))
      if(complex_ro) allocate(tmpi(maxstar))
      do istar=1,nstar
        mm=indstar(istar)-indstar(istar-1)
        tmp=0.d0
        if(complex_ro) tmpi=0.d0
        do iv=indstar(istar-1)+1,indstar(istar)
          do ig=1,ngroup
            call rotate(gbs(1,iv),gbs(2,iv),gbs(3,iv),x(1),x(2),x(3),
     &		            u(2,ig),2)
            do jv=indstar(istar-1)+1,indstar(istar)
              delt=abs(gbs(1,jv)-x(1))+abs(gbs(2,jv)-x(2))+
     &                                 abs(gbs(3,jv)-x(3))
              if(delt.lt.1.d-5) then
                jj=jv-indstar(istar-1)
                fas=-2.d0*pi*dot_product(x,shift(:,ig))
                cosn=cos(fas)
                sinn=sin(fas)
                tmp(jj)=tmp(jj)+cosn*real(fun(iv))
                if(complex_ro) then
                  tmp(jj)=tmp(jj)-sinn*imag(fun(iv))
                  tmpi(jj)=tmpi(jj)+cosn*imag(fun(iv))
     &                             +sinn*real(fun(iv))
                endif
              endif
            enddo  !!! over jv
          enddo  !!! over ig
        enddo  !!! over iv
        do iv=indstar(istar-1)+1,indstar(istar)
          ii=iv-indstar(istar-1)
          fun(iv)=tmp(ii)/ngroup
          if(complex_ro) fun(iv)=fun(iv)+(0.d0,1.d0)*tmpi(ii)/ngroup
        enddo  !!! over iv
      enddo  !!! over istar
      deallocate(tmp)
      if(complex_ro) deallocate(tmpi)
      end