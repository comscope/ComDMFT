      subroutine vgensp(rmax,nmax,a,a0,rout,nout,vect,nco,rco,lastv,
     &	                inf,err)
*****************************************************************
*       subroutine generates all vectors contained              *
*       in sphere rmax but not larger than nmax.                *
*                                                               *
*       input:                                                  *
*               rmax,nmax  - limiting parameters                *
*               a          - translations                       *
*               a0         - origin                             *
*       output:                                                 *
*               rout,      - sphere radii and number of all     *
*               nout         vectors in the sphere              *
*               vect       - vectors                            *
*               nco,rco(n) - number of co-sphers and their      *
*                            radii !r=rco(nco)                  *
*               lastv(n)   - number of last vector on           *
*                            co-sphere number n                 *
*               inf        - shows generation way: 0 - rmax,    *
*                            1 - nmax                           *
*****************************************************************
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: nmax
	real*8, intent(in) :: a(3,3),a0(3),rmax
	logical, intent(out) :: err
	integer, intent(out) :: nout,nco,lastv(nmax/2+1),inf
	real*8, intent(out) :: rout,rco(nmax/2),vect(3,nmax)
      integer :: nmax1,l1,l2,l3,nc,i,j,k,ivol,nc1m,nc1,ic,nv,last
      real*8 :: vs,s13,rrab,aa,bb,cc,vol,vol0,vv
	integer, allocatable :: iref(:),nvec(:),nco1(:),lastv1(:,:)
	real*8, allocatable :: vec(:,:),rw(:),rr1(:),rw1(:),vec1(:,:),
     &	                   rco1(:,:),vect1(:,:,:)
      err=.false.
      nmax1=max(nmax/numvol*4,10)
	allocate(vec(3,nmax))
	allocate(rw(nmax))
	allocate(iref(nmax))
	allocate(rr1(0:numvol))
	allocate(rw1(nmax1))
	allocate(vec1(3,nmax1))
	allocate(vect1(3,nmax1,numvol))
	allocate(rco1(nmax1,numvol))
	allocate(nvec(numvol))
	allocate(lastv1(nmax1,numvol))
	allocate(nco1(numvol))
      nout = 0
      rout = 0.0d0
c volume per site:
      vs=a(1,1)*(a(2,2)*a(3,3)-a(3,2)*a(2,3))
     &  +a(2,1)*(a(3,2)*a(1,3)-a(1,2)*a(3,3))
     &  +a(3,1)*(a(1,2)*a(2,3)-a(2,2)*a(1,3))
      vs=abs(vs)
c select rrab:
      s13=1.d0/3.d0
      rrab=(3.d0/4.d0/pi*vs*nmax)**s13
      rrab=min(rrab,rmax)
      inf=1
      if(rrab.eq.rmax) inf=0
c generation of all vectors in rrab sphere:
      aa=sqrt(a(1,1)*a(1,1)+a(2,1)*a(2,1)+a(3,1)*a(3,1))
      bb=sqrt(a(1,2)*a(1,2)+a(2,2)*a(2,2)+a(3,2)*a(3,2))
      cc=sqrt(a(1,3)*a(1,3)+a(2,3)*a(2,3)+a(3,3)*a(3,3))
      l1=rrab/vs*bb*cc+1
      l2=rrab/vs*aa*cc+1
      l3=rrab/vs*aa*bb+1
      nc=0
      do i=-l1,l1
         do j=-l2,l2
            do k=-l3,l3
               nc=nc+1
               if(nc.gt.nmax)then
                 err=.true.
                 goto 1
               endif
               vec(:,nc)=i*a(:,1)+j*a(:,2)+k*a(:,3)+a0
               rw(nc)=sqrt(vec(1,nc)**2+vec(2,nc)**2+vec(3,nc)**2)
               if((rw(nc)-rrab).ge.1.d-8) nc=nc-1
            enddo
         enddo
      enddo
c     *************************************************************
      if(nc.eq.1) then
         nco=1
         nout=0
         rout=0.d0
         rco(1)=0.d0
         lastv(1)=1
         vect(1,1)=0.d0
         vect(2,1)=0.d0
         vect(3,1)=0.d0
         goto 1
      endif
c     *************************************************************
      vol=4.d0*pi*rrab**3/3.d0
      vol0=vol/numvol
      vv=0.d0
      rr1(0)=-1.d0
      do ivol=1,numvol
         vv=vv+vol0
         rr1(ivol)=(3.d0*vv/4.d0/pi)**s13
      enddo
c     *************************************************************
      nc1m=0
      do ivol=1,numvol
         nc1=0
         do ic=1,nc
            if(rw(ic).gt.rr1(ivol-1).and.rw(ic).le.rr1(ivol)) then
               nc1=nc1+1
               rw1(nc1)=rw(ic)
               vec1(1,nc1)=vec(1,ic)
               vec1(2,nc1)=vec(2,ic)
               vec1(3,nc1)=vec(3,ic)
            endif
         enddo
         nc1m=max(nc1m,nc1)
c ordering:
         if(nc1.ne.0) call xorder(nc1,rw1,nco1(ivol),rco1(1,ivol),iref)
         if(nc1.eq.0) nco1(ivol)=0
c fulling of co-sphere:
c
         nvec(ivol)=nc1
         nv=0
         do i=1,nco1(ivol)   !!! over co-spheres for the sph layer IVOL
            do j=1,nc1       !!! over all vectors in the sph layer IVOL
               if(iref(j).eq.i) then  !! if the vector J in co-sphere I
                  nv=nv+1
                  do k=1,3
                     vect1(k,nv,ivol)=vec1(k,j)
                  enddo
               endif
            enddo
            lastv1(i,ivol)=nv
         enddo
      enddo
      if(nc1m.gt.nmax1)then
c         if(maswrk) write(iun,*)'<vgensph>: increase nmax '
c         call ending
        err=.true.
        goto 1
      endif
c output:
      nv=0
      nco=0
      last=0
      do ivol=1,numvol
         call dcopy(3*nvec(ivol),vect1(1,1,ivol),1,vect(1,nv+1),1)
         nv=nv+nvec(ivol)
         call dcopy(nco1(ivol),rco1(1,ivol),1,rco(nco+1),1)
         do i=1,nco1(ivol)
            nco=nco+1
            lastv(nco)=last+lastv1(i,ivol)
         enddo
         last=lastv(nco)
      enddo
      nout=nv
      rout=rco(nco)
1	deallocate(vec,rw,iref,rr1,rw1,vec1,vect1,rco1,nvec,lastv1,nco1)
      end
      
      
      subroutine vgensp2(rmax,nmax,a,b,rout,nout,vect,nco,rco,lastv,
     &	                 inf,err)
*****************************************************************
*       subroutine generates all vectors contained              *
*       in sphere rmax but not larger than nmax.                *
*                                                               *
*       input:                                                  *
*               rmax,nmax  - limiting parameters                *
*               a,b        - translations                       *
*       output:                                                 *
*               rout,      - sphere radii and number of all     *
*               nout         vectors in the sphere              *
*               vect       - vectors                            *
*               nco,rco(n) - number of co-sphers and their      *
*                            radii !r=rco(nco)                  *
*               lastv(n)   - number of last vector on           *
*                            co-sphere number n                 *
*               inf        - shows generation way: 0 - rmax,    *
*                            1 - nmax                           *
*****************************************************************
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: nmax
	real*8, intent(in) :: a(2),b(2),rmax
	logical, intent(out) :: err
	integer, intent(out) :: nout,nco,lastv(nmax/2+1),inf
	real*8, intent(out) :: rout,rco(nmax/2),vect(2,nmax)
      integer :: nmax1,l1,l2,nc,i,j,k,ivol,nc1m,nc1,ic,nv,last
      real*8 :: vs,rrab,aa,bb,vol,vol0,vv
	integer, allocatable :: iref(:),nvec(:),nco1(:),lastv1(:,:)
	real*8, allocatable :: vec(:,:),rw(:),rr1(:),rw1(:),vec1(:,:),
     &	                   rco1(:,:),vect1(:,:,:)
      err=.false.
      nmax1=max(nmax/numvol*4,10)
	allocate(vec(2,nmax))
	allocate(rw(nmax))
	allocate(iref(nmax))
	allocate(rr1(0:numvol))
	allocate(rw1(nmax1))
	allocate(vec1(2,nmax1))
	allocate(vect1(2,nmax1,numvol))
	allocate(rco1(nmax1,numvol))
	allocate(nvec(numvol))
	allocate(lastv1(nmax1,numvol))
	allocate(nco1(numvol))
c volume per site:
      vs=a(1)*b(2)-a(2)*b(1)
      vs=abs(vs)
c select rrab:
      rrab=sqrt(1.d0/pi*vs*nmax)
      rrab=min(rrab,rmax)
      inf=1
      if(rrab.eq.rmax) inf=0
c generation of all vectors in rrab sphere:
      aa=sqrt(a(1)*a(1)+a(2)*a(2))
      bb=sqrt(b(1)*b(1)+b(2)*b(2))
      l1=rrab/vs*bb+1
      l2=rrab/vs*aa+1
c
      nc=0
      do i=-l1,l1
        do j=-l2,l2
          nc=nc+1
          if(nc.gt.nmax)then
            err=.true.
            goto 1
          endif
          vec(1,nc)=i*a(1)+j*b(1)
          vec(2,nc)=i*a(2)+j*b(2)
          rw(nc)=sqrt(vec(1,nc)**2+vec(2,nc)**2)
          if((rw(nc)-rrab).ge.1.d-8) nc=nc-1
        enddo
      enddo
c     *************************************************************
      if(nc.eq.1) then
         nco=1
         nout=0
         rout=0.d0
         rco(1)=0.d0
         lastv(1)=1
         vect(1,1)=0.d0
         vect(2,1)=0.d0
         goto 1
      endif
c     *************************************************************
      vol=pi*rrab**2
      vol0=vol/numvol
      vv=0.d0
      rr1(0)=-1.d0
      do ivol=1,numvol
         vv=vv+vol0
         rr1(ivol)=sqrt(vv/pi)
      enddo
c     *************************************************************
      nc1m=0
      do ivol=1,numvol
         nc1=0
         do ic=1,nc
            if(rw(ic).gt.rr1(ivol-1).and.rw(ic).le.rr1(ivol)) then
               nc1=nc1+1
               rw1(nc1)=rw(ic)
               vec1(1,nc1)=vec(1,ic)
               vec1(2,nc1)=vec(2,ic)
            endif
         enddo
         nc1m=max(nc1m,nc1)
c ordering:
         if(nc1.ne.0) call xorder(nc1,rw1,nco1(ivol),rco1(1,ivol),iref)
         if(nc1.eq.0) nco1(ivol)=0
c fulling of co-sphere:
c
         nvec(ivol)=nc1
         nv=0
         do i=1,nco1(ivol)   !!! over co-spheres for the sph layer IVOL
            do j=1,nc1       !!! over all vectors in the sph layer IVOL
               if(iref(j).eq.i) then  !! if the vector J in co-sphere I
                  nv=nv+1
                  do k=1,2
                     vect1(k,nv,ivol)=vec1(k,j)
                  enddo
               endif
            enddo
            lastv1(i,ivol)=nv
         enddo
      enddo
      if(nc1m.gt.nmax1)then
c         if(maswrk) write(iun,*)'<vgensph>: increase nmax '
c         call ending
        err=.true.
        goto 1
      endif
c output:
      nv=0
      nco=0
      last=0
      do ivol=1,numvol
         call dcopy(2*nvec(ivol),vect1(1,1,ivol),1,vect(1,nv+1),1)
         nv=nv+nvec(ivol)
         call dcopy(nco1(ivol),rco1(1,ivol),1,rco(nco+1),1)
         do i=1,nco1(ivol)
            nco=nco+1
            lastv(nco)=last+lastv1(i,ivol)
         enddo
         last=lastv(nco)
      enddo
      nout=nv
      rout=rco(nco)
1	deallocate(vec,rw,iref,rr1,rw1,vec1,vect1,rco1,nvec,lastv1,nco1)
      end
