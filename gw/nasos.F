      subroutine nasos(wr,rmt)
*********************************************************
*       subroutine for determining mt-sphere            *
*       radii according to  input weight.               *
*********************************************************
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	real*8, intent(in) :: wr(natom)
	real*8, intent(out) :: rmt(natom)
      integer :: nvcmax,i,ns,j,k,iat,ncheck
	real*8 :: rw,dog,al0,dij,al,g
	integer, allocatable :: flag(:),isat(:)
	real*8, allocatable :: tauc(:,:),site(:,:)
	nvcmax=30000
	allocate(flag(natom))
	allocate(isat(nvcmax))
	allocate(tauc(3,natom))
	allocate(site(3,nvcmax))
c generation all sites in 8 cell:
      do i=1,natom
        call reduc(rbas(1,1),rbas(1,2),rbas(1,3),tau(1,i),tauc(1,i),1)
      enddo
      ns=0
      do i=0,1
         do j=0,1
            do k=0,1
               do iat=1,natom
                  ns=ns+1
                  if(ns.gt.nvcmax)then
                     if(maswrk) write(iun,*) '<NASOS>: increase nvcmax'
                     call ending
                  endif
                  isat(ns)=iat
	            site(:,ns)=tauc(:,iat)-i*rbas(:,1)-j*rbas(:,2)
     &				                              -k*rbas(:,3)
               enddo
            enddo
         enddo
      enddo
c basis set:
      do i=1,natom
        call reduc(rbas(1,1),rbas(1,2),rbas(1,3),tau(1,i),tauc(1,i),0)
      enddo
      do iat=1,natom
         do i=1,ns
            rw=(site(1,i)-tauc(1,iat))**2+
     &         (site(2,i)-tauc(2,iat))**2+
     &         (site(3,i)-tauc(3,iat))**2
            if(rw.lt.1.d-08) then
               do k=1,3
                  site(k,i)=site(k,iat)
               enddo
               do k=1,3
                  site(k,iat)=tauc(k,iat)
               enddo
               goto 110
            endif
         enddo
110      continue
      enddo
c start up:
      do i=1,natom
        rmt(i)=0.d0
        flag(i)=0
      enddo
c blowing coefficient:
      dog=0
125   al0=1.d10
      dog=dog+1
      if (dog.gt.100) then
         if(maswrk) write(iun,*)'stop from nasos: implicit cycle abort'
         call ending
      endif
      do i=1,natom
         if(flag(i).eq.1) goto 140
         do j=1,ns
            if(i.eq.j) goto 130
            dij=sqrt((site(1,i)-site(1,j))**2+(site(2,i)-site(2,j))**2+
     &               (site(3,i)-site(3,j))**2)
            al=(dij-rmt(isat(j))*flag(isat(j)))/
     &         (wr(i)+wr(isat(j))*(1-flag(isat(j))))
            al0=min(al0,al)
130         continue
         enddo
140      continue
      enddo
c blowing:
      do i=1,natom
         if(flag(i).eq.0) rmt(i)=wr(i)*al0
      enddo
C flag set:
      do i=1,natom
         do j=1,ns
            if(i.eq.j) goto 160
            dij=sqrt((site(1,i)-site(1,j))**2+(site(2,i)-site(2,j))**2+
     &               (site(3,i)-site(3,j))**2)
            g=abs(dij-rmt(i)-rmt(isat(j)))
            if(g.le.1.d-08) then
               flag(i)=1
               goto 170
            endif
160         continue
         enddo
170      continue
      enddo
c check:
      ncheck=0
      do i=1,natom
         ncheck=ncheck+flag(i)
      enddo
      if(ncheck.lt.natom) goto 125
	deallocate(flag,isat,tauc,site)
      end