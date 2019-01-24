      subroutine poiss1(vc1,vc2,l,nrad,h,r,weight,rl,rotot)
*************************************************************
* The solution of the Poisson's equation within mt-spheres. *
*************************************************************
	use manager_mod
      implicit none
	integer, intent(in) :: l,nrad
      real*8, intent(in) :: rotot(0:nrad),r(0:nrad),weight(0:nrad),
     &                      rl(0:nrad),h
      real*8, intent(out) :: vc1(0:nrad),vc2(0:nrad)
	integer :: irad,l1
	real*8, allocatable :: work1(:),work2(:)
c   initiation
	allocate(work2(0:nrad),work1(0:nrad))
c   solution of the poisson equation
      do irad=0,nrad
        work1(irad)=rotot(irad)*weight(irad)*rl(irad)
      enddo
      call intalt(nrad+1,work1*h,work2)
      vc1(0)=0.0
      do irad=1,nrad
        vc1(irad)=work2(irad)
        do l1=1,l+1
          vc1(irad)=vc1(irad)/r(irad)
	  enddo
      enddo
      work1(0)=0.0
      do irad=1,nrad
        work1(irad)=rotot(irad)*weight(irad)
        do l1=1,l+1
          work1(irad)=work1(irad)/r(irad)
	  enddo
      enddo
      call intalt(nrad+1,work1*h,work2)
      vc2(0)=0.0
      do irad=1,nrad
        vc2(irad)=(work2(nrad)-work2(irad))*rl(irad)
	enddo
      if(l.eq.0)vc2(0)=work2(nrad)
	deallocate(work2,work1)
      end
