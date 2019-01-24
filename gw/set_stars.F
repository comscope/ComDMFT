      subroutine set_stars(nco,nsta,indsta,qabs_sta,ii,err)
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: nco,ii
      logical, intent(out) :: err
      integer, intent(out) :: nsta,indsta(0:ii*nco)
      real*8, intent(out) :: qabs_sta(ii*nco)
      integer :: ico,i,ig,j,l,k
      real*8 :: x(3),delt,gb
c     &&&&&&&&&&   Generation of Stars &&&&&&&&&&&&&&&&&&&&&&&&
c     indstar(istar) - last index star=istar
      nsta=0
      indsta(0)=0
      err=.false.
      do ico=1,nco
        i=lastg(ico-1)+1
8       nsta=nsta+1
        if(nsta>ii*nco) then
          err=.true.
          return
        endif
        qabs_sta(nsta)=rcosph(ico)
        indsta(nsta)=i
        do ig=1,ngroup
          call rotate(gbs(1,i),gbs(2,i),gbs(3,i),x(1),x(2),x(3),
     &                u(2,ig),2)
          do j=indsta(nsta-1)+1,indsta(nsta)
            delt=abs(gbs(1,j)-x(1))+abs(gbs(2,j)-x(2))+
     &                              abs(gbs(3,j)-x(3))
            if(delt.lt.1.d-5) goto 7
          enddo
          indsta(nsta)=indsta(nsta)+1
          do j=indsta(nsta),lastg(ico)
            delt=abs(gbs(1,j)-x(1))+abs(gbs(2,j)-x(2))+
     &                              abs(gbs(3,j)-x(3))
            if(delt.lt.1.d-5) then
              l=indsta(nsta)
              if(j.ne.l) then
                do k=1,3
                  gb=gbs(k,l)
                  gbs(k,l)=gbs(k,j)
                  gbs(k,j)=gb
                enddo
              endif
              goto 7
            endif
          enddo
7         continue
        enddo   !!! over ig
        if(indsta(nsta).lt.lastg(ico)) then
          i=indsta(nsta)+1
          goto 8
        endif
      enddo  !!! over ico
      end
      
      
      subroutine set_stars2d(nco,nsta,indsta,qabs_sta,ii,err)
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: nco,ii
      logical, intent(out) :: err
      integer, intent(out) :: nsta,indsta(0:ii*nco)
      real*8, intent(out) :: qabs_sta(ii*nco)
      integer :: ico,i,ig,j,l,k
      real*8 :: x(3),delt,gb
c     &&&&&&&&&&   Generation of Stars &&&&&&&&&&&&&&&&&&&&&&&&
c     indstar(istar) - last index star=istar
      nsta=0
      indsta(0)=0
      err=.false.
      do ico=1,nco
        i=lastg(ico-1)+1
8       nsta=nsta+1
        if(nsta>ii*nco) then
          err=.true.
          return
        endif
        qabs_sta(nsta)=rcosph(ico)
        indsta(nsta)=i
        do ig=1,ngroup
          call rotate2d(gbs(1,i),gbs(2,i),x(1),x(2),u_opt(1,1,ig),2)
          do j=indsta(nsta-1)+1,indsta(nsta)
            delt=abs(gbs(1,j)-x(1))+abs(gbs(2,j)-x(2))
            if(delt.lt.1.d-5) goto 7
          enddo
          indsta(nsta)=indsta(nsta)+1
          do j=indsta(nsta),lastg(ico)
            delt=abs(gbs(1,j)-x(1))+abs(gbs(2,j)-x(2))
            if(delt.lt.1.d-5) then
              l=indsta(nsta)
              if(j.ne.l) then
                do k=1,2
                  gb=gbs(k,l)
                  gbs(k,l)=gbs(k,j)
                  gbs(k,j)=gb
                enddo
              endif
              goto 7
            endif
          enddo
7         continue
        enddo   !!! over ig
        if(indsta(nsta).lt.lastg(ico)) then
          i=indsta(nsta)+1
          goto 8
        endif
      enddo  !!! over ico
      end