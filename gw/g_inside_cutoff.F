      subroutine g_inside_cutoff(num,pt,cut,order,ind,indx,n,n0)
c     num - max index of GBS vector to look ---
c     pt(3) - k-point
c     cut   - cutoff
c     order=T if we order the vectors (ascending)
c     ind=T if we need index array INDX(n0)      
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      logical, intent(in) :: order,ind
      integer, intent(in) :: num,n0
      real*8, intent(in) :: pt(3),cut
      integer, intent(out) :: n,indx(n0)
	integer :: ig,i0,i1,i
	real*8 :: ek,gk,p(3),gkmax
	integer, allocatable :: flag(:)
	real*8, allocatable :: gb(:,:)
	allocate(gb(2,n0))
	allocate(flag(num))
	n=0
c  ---- First include all vectors strictly inside CUTOFF --------
      gkmax=0.d0
      flag=0
      do ig=1,num
        p=pt+gbs(:,ig)
        ek=dot_product(p,p)
        gk=sqrt(ek)
        if(gk<=cut) then
          gkmax=max(gk,gkmax)
          flag(ig)=1
          n=n+1
          if(order) then
            do i0=1,n-1
              if(gb(2,i0)>gk) then
                do i1=n,i0+1,-1
                  do i=1,2
                    gb(i,i1)=gb(i,i1-1)
                  enddo
                enddo
                gb(1,i0)=dfloat(ig)
                gb(2,i0)=gk
                goto 1
              endif
            enddo
          endif
          gb(1,n)=dfloat(ig)
          gb(2,n)=gk
        endif
1       continue
      enddo
c  ---- Now we add vectors which are very close to the max included ---
      do ig=1,num
        if(flag(ig)==1) cycle
        p=pt+gbs(:,ig)
        ek=dot_product(p,p)
        gk=sqrt(ek)
        if(gk-gkmax<1.d-10) then
          n=n+1
          if(order) then
            do i0=1,n-1
              if(gb(2,i0)>gk) then
                do i1=n,i0+1,-1
                  do i=1,2
                    gb(i,i1)=gb(i,i1-1)
                  enddo
                enddo
                gb(1,i0)=dfloat(ig)
                gb(2,i0)=gk
                goto 2
              endif
            enddo
          endif
          gb(1,n)=dfloat(ig)
          gb(2,n)=gk
        endif
2       continue
      enddo
c ------------------------------------------------------
      if(ind) indx(1:n)=idnint(gb(1,1:n))
	deallocate(gb,flag)
      end