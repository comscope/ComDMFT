      subroutine sitegen(rmax,a,b,c,ta,nat,ns,site,iref,nvcmax)
*****************************************************************
*       subroutine for generation of all sites contained        *
*       in sphere rmax.                                         *
*                                                               *
*       input:                                                  *
*               a,b,c      - translations                       *
*               tau,natom  - basis                              *
*       output:                                                 *
*               ns         - number of all sites in the sphere  *
*               site       - vectors                            *
*               iref(n)    - shows the site sort                *
*****************************************************************
	use manager_mod
	use parallel_mod
      use sitegen_mod
      implicit none
      integer :: nvcmax,iref(nvcmax),nat,ns,i,nvec,nco,inf,j,k
      real*8 :: a(3),b(3),c(3),site(3,nvcmax),ta(3,nat),radd,s,rrab,
     &	        rmax,rou,sch
      radd=0
      do i=1,nat
        call reduc(a,b,c,ta(1,i),tauc_site(1,i),0)
        s=sqrt(tauc_site(1,i)**2+tauc_site(2,i)**2+tauc_site(3,i)**2)
        radd=max(radd,s)
	enddo
c select rrab:
      rrab=rmax+radd
c generation of all vectors in rrab sphere:
      call vgensph(rrab,nvcmax,a,b,c,rou,nvec,vect_site,nco,rco_site,
     &             lastg_site,inf,nvcmax)
c     ********************************************************
      if(inf.eq.1.and.maswrk) write(iun,*)'warning
     &  from sitegen: some sites may be missed'
c
      ns=0
      do i=0,nvec
        do j=1,nat
          ns=ns+1
          if(ns.gt.nvcmax)then
            if(maswrk) write(iun,*)'<sitegen>: ns>nvcmax'
            call ending
          endif
          iref(ns)=j
          do k=1,3
            site(k,ns)=vect_site(k,i)+tauc_site(k,j)
          enddo
          sch=sqrt(site(1,ns)**2+site(2,ns)**2+site(3,ns)**2)
          if((sch-rmax).gt.1.d-08) ns=ns-1
        enddo
      enddo
      end