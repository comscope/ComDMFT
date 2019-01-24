      subroutine vgensph(rmax,nmax,a,b,c,r,n,vect,nco,rco,lastg,inf,
     &                   nvcmax)
*****************************************************************
*       subroutine generates all vectors contained              *
*       in sphere rmax but not larger than nmax.                *
*                                                               *
*       input:                                                  *
*               rmax,nmax  - limiting parameters                *
*               a,b,c      - translations                       *
*       output:                                                 *
*               r,n        - sphere radii and number of all     *
*                            vectors in the sphere              *
*               vect       - vectors
*               nco,rco(n) - number of co-sphers and there      *
*                            radii !r=rco(nco)                  *
*               lastg(n)   - number of last vector on           *
*                            co-sphere number n                 *
*               inf        - shows generation way: 0 - rmax,    *
*                            1 - nmax                           *
*****************************************************************
	use manager_mod
      use parallel_mod
        use sitegen_mod
	use units_mod
      implicit none
      integer :: nvcmax,lastg(nvcmax/2),nmax,inf,l1,l2,l3,nc,i,j,k,nco,
     &	       n,nvec
      real*8 :: a(3),b(3),c(3),rco(nvcmax/2),vect(3,0:nvcmax),vs,rrab,
     &	      rmax,drrab,dog,aa,bb,cc,r
c volume per site:
        vs=a(1)*(b(2)*c(3)-b(3)*c(2))+a(2)*(b(3)*c(1)-b(1)*c(3))+
     &    a(3)*(b(1)*c(2)-b(2)*c(1))
        vs=abs(vs)
c select rrab:
        rrab=(3.d0/4.d0/pi*vs*nmax)**(1./3.)
        rrab=min(rrab,rmax)
        inf=1
        if(rrab.eq.rmax) inf=0
        drrab=0
        dog=0
c generation of all vectors in rrab sphere:
 100    rrab=rrab+drrab
        if (dog.gt.100) then
           if(maswrk) write(iun,*)'vgensph: implicit cycle abort'
           call ending
        endif
        aa=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
        bb=sqrt(b(1)*b(1)+b(2)*b(2)+b(3)*b(3))
        cc=sqrt(c(1)*c(1)+c(2)*c(2)+c(3)*c(3))
        l1=rrab/vs*bb*cc+1
        l2=rrab/vs*aa*cc+1
        l3=rrab/vs*aa*bb+1
c
                nc=0
                do 30 i=-l1,l1
                do 20 j=-l2,l2
                do 10 k=-l3,l3
                nc=nc+1
                if(nc.gt.nvcmax)then
                  if(maswrk) write(iun,*)'<vgensph>: increase nvcmax '
                  call ending
                endif
                vec_vgen(1,nc)=i*a(1)+j*b(1)+k*c(1)
                vec_vgen(2,nc)=i*a(2)+j*b(2)+k*c(2)
                vec_vgen(3,nc)=i*a(3)+j*b(3)+k*c(3)
                rw_vgen(nc)=sqrt(vec_vgen(1,nc)**2+vec_vgen(2,nc)**2
     &                     +vec_vgen(3,nc)**2)
                if((rw_vgen(nc)-rrab).ge.1.d-08) nc=nc-1
  10            continue
  20            continue
  30            continue
        if(nc.gt.nvcmax) then
        drrab=-rrab*(nc-nvcmax)/3.d0/nvcmax
        dog=dog+1
        goto 100
        endif
c inserted 12.11.92
        if(nc.eq.1) then
        nco=1
        n=0
        r=0.
        rco(1)=0.d0
        lastg(1)=0
        vect(1,0)=0.d0
        vect(2,0)=0.d0
        vect(3,0)=0.d0
        return
        endif
c ordering:
        call xorder(nc,rw_vgen,nco,rco,iref_vgen)
c fulling of co-sphere:
c
        nvec=0
        do 60 i=2,nco
        rco(i-1)=rco(i)
        do 50 j=1,nc
        if(iref_vgen(j).eq.i) then
        nvec=nvec+1
        do 40 k=1,3
        vect(k,nvec)=vec_vgen(k,j)
  40    continue
        endif
  50    continue
        lastg(i-1)=nvec
        if(nvec.gt.nmax) goto 70
  60    continue
c output:
        nco=nco-1
        n=nvec
        goto 80
  70    nco=i-2
        n=lastg(i-2)
  80    r=rco(nco)
      vect(1,0)=0.d0
      vect(2,0)=0.d0
      vect(3,0)=0.d0
      end
