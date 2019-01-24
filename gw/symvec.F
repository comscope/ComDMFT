      subroutine symvec(fun,jreal,jimag)
c     nplwro  - number of plane wave for density (in sphere cutoffro)
	use solid_mod
	use units_mod
      implicit none
      logical, intent(in) :: jreal,jimag
      complex*16, intent(out) :: fun(3,nplwro)
	integer :: istar,mm,iv,ig,jv,jj,i,j,ii
      real*8 :: v1(3,2),vv(3),delt,fas,cosn,sinn,ddot
	real*8, allocatable :: tmp(:,:),tmpi(:,:)
      if(inv<2) allocate(tmp(3,maxstar))
      if(inv/=1) allocate(tmpi(3,maxstar))
      do istar=1,nstar
        mm=3*(indstar(istar)-indstar(istar-1))
        if(jreal) call dcopy(mm,0.d0,0,tmp,1)
        if(jimag) call dcopy(mm,0.d0,0,tmpi,1)
        do iv=indstar(istar-1)+1,indstar(istar)
          do ig=1,ngroup
            call rotate(gbs(1,iv),gbs(2,iv),gbs(3,iv),vv(1),
     &                  vv(2),vv(3),u(2,ig),2)
            do jv=indstar(istar-1)+1,indstar(istar)
              delt=abs(gbs(1,jv)-vv(1))+abs(gbs(2,jv)-vv(2))+
     &                                     abs(gbs(3,jv)-vv(3))
              if(delt.lt.1.d-5) then
                jj=jv-indstar(istar-1)
                fas=-2.d0*pi*ddot(3,vv,1,shift(1,ig),1)
                cosn=cos(fas)*ifia(ig)
                sinn=sin(fas)*ifia(ig)
                do i=1,3
                  if(jreal) v1(i,1)=0.d0
                  if(jimag) v1(i,2)=0.d0
                  do j=1,3
                    if(jreal) v1(i,1)=v1(i,1)+agm(i,j,ig)
     &                                       *real(fun(j,iv))
                    if(jimag) v1(i,2)=v1(i,2)+agm(i,j,ig)
     &                                       *imag(fun(j,iv))
                  enddo
                  if(jreal) tmp(i,jj)=tmp(i,jj)+cosn*v1(i,1)
     &				                         -sinn*v1(i,2)
                  if(jimag) tmpi(i,jj)=tmpi(i,jj)+cosn*v1(i,2)
     &				                           +sinn*v1(i,1)
                enddo
              endif
            enddo  !!! over jv
          enddo  !!! over ig
        enddo  !!! over iv
        do iv=indstar(istar-1)+1,indstar(istar)
          ii=iv-indstar(istar-1)
          do i=1,3
            if(jreal) fun(i,iv)=tmp(i,ii)/ngroup
            if(jimag) fun(i,iv)=(0.d0,1.d0)*tmpi(i,ii)/ngroup
          enddo
        enddo  !!! over iv
      enddo  !!! over istar
      if(inv<2) deallocate(tmp)
      if(inv/=1) deallocate(tmpi)
      end