      subroutine parsbl(t,nt,i)
      implicit none
      integer :: nt,i
      character*1 t(nt)
  99  if(t(i).ne.'_') return
      i=i+1
      if(i.gt.nt) return
      goto 99
      end
      
      subroutine parsgn(t,nt,g,g_j,a,ifiaa,lmax,jmax,ng,h,h_j,su1,su1_j,
     &                  fac)
c  Parses the string t containing the group element describers
c  and makes a matrix g and vector a for each element.
c  space group ops are defined as:  (g,a)(p)=g*p+a
c
c  Point group element describers:
c    Rnv = n-fold rot around vec v, n pos integer   ) also products
c    Mv  = mirror taking vec v into -v              ) of form 'a*b' of
c    I   = inversion                                ) any two
c    -   = J => -J                                  ) of these
c
c  Glide planes, screw axes etc:
c    append  ':Tv'  to any point operation above where v is the
c    translation vector for the operation.
c
c  Vectors v are all given as (x,y,z) with real numbers x,y,z
c  or abbreviate   'D' for (1,1,1), 'X' for (1,0,0),  Y,Z  similar.
c
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer ifiaa(*),nt,i,ng,lmax,jmax,l,iwig,ist,ii,jj
      character*1 t(nt)
      real*8 g(maxwig,*),h(maxwig),a(3,*),su1(-lmax:lmax,-lmax:lmax),
     &       fac(0:*)
      complex*16 g_j(maxwig_j,*),su1_j(jmax+1,jmax+1),
     &           h_j(maxwig_j)
      ng=0
      i=1
  90  call parsbl(t,nt,i)
      if(i.gt.nt) return
      ng=ng+1
      ifiaa(ng)=1
      if(t(i).eq.'-') then
        ifiaa(ng)=-1
        i=i+1
      endif
      call parsop(t,nt,i,g(1,ng),g_j(1,ng),lmax,jmax,su1,su1_j,fac)
      if(t(i).eq.'*') then
        i=i+1
        call parsop(t,nt,i,h,h_j,lmax,jmax,su1,su1_j,fac)
        do l=0,lmax
          iwig=l*(2*l-1)*(2*l+1)/3+1
          call dcopy((2*l+1)**2,g(iwig,ng),1,su1(-lmax,-lmax),1)
          call prd(su1(-lmax,-lmax),h(iwig),g(iwig,ng),l)
        enddo
        if(irel.eq.2) then
          do l=0,lmax
            ist=-1
            if(l.eq.0) ist=1
            do ii=ist,1,2
              jj=l+l+ii
              iwig=4*l*(2*l*l+1)/3+2*(ii+1)*l*l+1
              call zcopy((jj+1)**2,g_j(iwig,ng),1,su1_j,1)
              call prd_j(su1_j,h_j(iwig),g_j(iwig,ng),jj)
            enddo
          enddo
        endif
      endif
      call dcopy(3,0.d0,0,a(1,ng),1)
      if(t(i).eq.':') then
        i=i+1
        if(t(i).ne.'T') then
          if(maswrk) write(iun,*)'*** parsgn: "T" expected after ":"'
          call ending
        endif
        i=i+1
        call parsvc(t,nt,i,a(1,ng),rbas)
      endif
      if(t(i).ne.'_') then
        if(maswrk) write(iun,*)'*** no blank found after op'
        call ending
      endif
      goto 90
      end
      
      subroutine parsop(t,nt,i,a,djgen,lmax,jmax,su1,su1_j,fac)
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer :: nt,i,k,l,m,m1,lmax,jmax,ndig,iwig,ist,ii,iwig0,jj,
     &           im,im1,ik
      integer*8 :: nrot
      character*1 t(nt)
      real*8 :: v(3),a(maxwig),su1(-lmax:lmax,-lmax:lmax),fac(0:*),sp,
     &          v12,teta,fi,om,ddot,su,ss,wign
      complex*16 :: wign_j,djgen(maxwig_j),su1_j(jmax+1,jmax+1),cu,
     &              suj
      if(t(i).eq.'R') then
        i=i+1
        call parspi(t,nt,i,nrot,ndig)
        if(ndig.eq.0) then
          if(maswrk) write(iun,*)'*** no integer found after "R"'
          call ending
        endif
        call parsvc(t,nt,i,v,rbas)
        sp=ddot(3,v,1,v,1)
        sp=1.d0/dsqrt(sp)
        do k=1,3
          v(k)=v(k)*sp
        enddo
        v12=v(1)**2+v(2)**2
        if(v12.gt.1.d-12) then
          teta=acos(v(3))
          fi=acos(v(1)/sqrt(v12))
          if(v(2).lt.0.d0) fi=-fi
        else
          teta=0.d0
          fi=0.d0
        endif
        om=2.d0*pi/nrot
        do l=0,lmax
          do m=-l,l
            do m1=-l,l
              su1(m,m1)=wign(l,m,m1,fi,teta,0.d0,fac)
            enddo
          enddo
          do m=-l,l
            do m1=-l,l
              iwig=l*(2*l-1)*(2*l+1)/3+(2*l+1)*(l+m1)+l+m+1
              su=0.d0
              do k=-l,l
                ss=k*om
                su=su+su1(m1,k)*(su1(m,k)*dcos(ss)+su1(m,-k)*dsin(ss))
              enddo
              a(iwig)=su
            enddo
          enddo
        enddo
        if(irel.eq.2) then
          do l=0,lmax
            ist=-1
            if(l.eq.0) ist=1
            do ii=ist,1,2
              iwig0=(8*l**3+12*l**2+10*l+3)/3+l*(2*l+1)*ii
              jj=l+l+ii
              do m=-jj,jj,2
                im=(m+jj)/2+1
                do m1=-jj,jj,2
                  im1=(m1+jj)/2+1
                  su1_j(im,im1)=wign_j(jj,m,m1,fi,teta,0.d0,fac)
                enddo
              enddo
              do m=-jj,jj,2
                im=(m+jj)/2+1
                do m1=-jj,jj,2
                  im1=(m1+jj)/2+1
                  iwig=iwig0+(l+l+ii+1)*(m1+ii)/2+(ii+m)/2
                  suj=(0.d0,0.d0)
                  do k=-jj,jj,2
                    ik=(k+jj)/2+1
                    ss=k*om/2.d0
                    cu=dcmplx(0.d0,ss)
                    cu=exp(cu)
                    suj=suj+dconjg(su1_j(im1,ik))*su1_j(im,ik)*cu
                  enddo
                  djgen(iwig)=suj
                enddo
              enddo
            enddo
          enddo
        endif
      else if(t(i).eq.'M') then
        i=i+1
        call parsvc(t,nt,i,v,rbas)
        sp=ddot(3,v,1,v,1)
        sp=1.d0/dsqrt(sp)
        do k=1,3
          v(k)=v(k)*sp
        enddo
        v12=v(1)**2+v(2)**2
        if(v12.gt.1.d-12) then
          teta=acos(v(3))
          fi=acos(v(1)/sqrt(v12))
          if(v(2).lt.0.d0) fi=-fi
        else
          teta=0.d0
          fi=0.d0
        endif
        do l=0,lmax
          do m=-l,l
            do m1=-l,l
              su1(m,m1)=wign(l,m,m1,fi,teta,0.d0,fac)
            enddo
          enddo
          do m=-l,l
            do m1=-l,l
              iwig=l*(2*l-1)*(2*l+1)/3+(2*l+1)*(l+m1)+l+m+1
              su=0.d0
              do k=-l,l
                su=su+(-1.d0)**(l+k)*su1(m1,k)*su1(m,k)
              enddo
              a(iwig)=su
            enddo
          enddo
        enddo
        if(irel.eq.2) then
          do l=0,lmax
            ist=-1
            if(l.eq.0) ist=1
            do ii=ist,1,2
              iwig0=(8*l**3+12*l**2+10*l+3)/3+l*(2*l+1)*ii
              jj=l+l+ii
              do m=-jj,jj,2
                im=(m+jj)/2+1
                do m1=-jj,jj,2
                  im1=(m1+jj)/2+1
                  su1_j(im,im1)=wign_j(jj,m,m1,fi,teta,0.d0,fac)
                enddo
              enddo
              do m=-jj,jj,2
                im=(m+jj)/2+1
                do m1=-jj,jj,2
                  im1=(m1+jj)/2+1
                  iwig=iwig0+(l+l+ii+1)*(m1+ii)/2+(ii+m)/2
                  suj=(0.d0,0.d0)
                  do k=-jj,jj,2
                    ik=(k+jj)/2+1
                    suj=suj+dconjg(su1_j(im1,ik))*su1_j(im,ik)
     &                 *(-1.d0)**(l+(k+1)/2)*(0.d0,1.d0)
                  enddo
                  djgen(iwig)=suj
                enddo
              enddo
            enddo
         enddo
        endif
      else if(t(i).eq.'I') then
        do l=0,lmax
          do m=-l,l
            do m1=-l,l
              iwig=l*(2*l-1)*(2*l+1)/3+(2*l+1)*(l+m1)+l+m+1
              if (m.eq.m1)then
                a(iwig)=(-1)**l
              else
                a(iwig)=0.d0
              endif
            enddo
          enddo
        enddo
        if(irel.eq.2) then
          do l=0,lmax
            ist=-1
            if(l.eq.0) ist=1
            do ii=ist,1,2
              iwig0=(8*l**3+12*l**2+10*l+3)/3+l*(2*l+1)*ii
              jj=l+l+ii
              do m=-jj,jj,2
                do m1=-jj,jj,2
                  iwig=iwig0+(l+l+ii+1)*(m1+ii)/2+(ii+m)/2
                  if (m.eq.m1)then
                    ss=(-1)**l
                  else
                    ss=0.d0
                  endif
                  djgen(iwig)=dcmplx(ss,0.d0)
                enddo
              enddo
            enddo
          enddo
        endif
        i=i+1
      else if(t(i).eq.'E') then
        do l=0,lmax
          do m=-l,l
            do m1=-l,l
              iwig=l*(2*l-1)*(2*l+1)/3+(2*l+1)*(l+m1)+l+m+1
              if (m.eq.m1)then
                a(iwig)=1.d0
              else
                a(iwig)=0.d0
              endif
            enddo
          enddo
        enddo
        if(irel.eq.2) then
          do l=0,lmax
            ist=-1
            if(l.eq.0) ist=1
            do ii=ist,1,2
              iwig0=(8*l**3+12*l**2+10*l+3)/3+l*(2*l+1)*ii
              jj=l+l+ii
              do m=-jj,jj,2
                do m1=-jj,jj,2
                  iwig=iwig0+(l+l+ii+1)*(m1+ii)/2+(ii+m)/2
                  if (m.eq.m1)then
                    ss=1.d0
                  else
                    ss=0.d0
                  endif
                  djgen(iwig)=dcmplx(ss,0.d0)
                enddo
              enddo
            enddo
          enddo
        endif
        i=i+1
      else
        if(maswrk) write(iun,*)'*** op must start with "R","M" or "I"'
        call ending
      endif
      end
      
      
      subroutine parspi(t,nt,i,integr,ndig)
      use manager_mod
      implicit none
      integer :: nt,i,ndig,is,mxdig
      integer*8 :: integr
      character*1 t(nt),csym(10)
      data csym/'0','1','2','3','4','5','6','7','8','9'/
      ! mxdig: the maximum number of digits that fit into integer*8
      mxdig=log10(dble(huge(integr)))
      integr=0
      ndig=0
      i=i-1
  99  i=i+1
      do 1 is=1,10
        if(csym(is).eq.t(i)) then
          integr=10*integr+is-1
          ndig=ndig+1
          goto 99
        endif
  1   continue
      if (ndig.gt.mxdig) then
        write(iun,*)'Subroutine parspi: Integer overflow detected'
        write(iun,*)'Too many digits:'
        write(iun,*)'text    =',t
        write(iun,*)'position=',i
        write(iun,*)'mxdigits=',mxdig
        write(iun,*)'# digits=',ndig
        write(iun,*)'integer =',integr
        call ending
      endif
      end
      
      subroutine parssr(t,nt,i,r)
      implicit none
      integer :: nt,i,isig,ndig
      integer*8 :: intg1, intg2
      character*1 t(nt)
      real*8 r
      isig=1
      if(t(i).eq.'-') isig=-1
      if(t(i).eq.'-'.or.t(i).eq.'+') i=i+1
      call parspi(t,nt,i,intg1,ndig)
      if(t(i).eq.'.') then
        i=i+1
        call parspi(t,nt,i,intg2,ndig)
        r=isig*(intg1+intg2*(0.1d0**ndig))
      else
        r=intg1*isig
      endif
      end
      
      subroutine parsvc(t,nt,i,v,rbas)
      use manager_mod
      use parallel_mod
      implicit none
      integer :: nt,i,ii,ndig,m,jj
      integer*8 :: integr
      ! if ibracket==1: vectors given in terms of lattice vectors
      ! if ibracket==2: vectors given in terms of cartesian coordinates
      integer :: ibracket
      real*8 :: r,v(3),rbas(3,3),fac
      character*1 t(nt)
      ibracket=0
      v(1)=0.d0
      v(2)=0.d0
      v(3)=0.d0
      if(t(i).eq.'X') v(1)=1.d0
      if(t(i).eq.'Y') v(2)=1.d0
      if(t(i).eq.'Z') v(3)=1.d0
      if(t(i).eq.'X'.or.t(i).eq.'Y'.or.t(i).eq.'Z') goto 90
      if(t(i).eq.'D') then
        v(1)=1.d0
        v(2)=1.d0
        v(3)=1.d0
        goto 90
      endif
      if(t(i).eq.'S') then
        v(1)=0.d0
        v(2)=0.d0
        v(3)=1.d0
        i=i+1
        call parspi(t,nt,i,integr,ndig)
        if(integr.ge.10) then
          fac=dfloat(mod(integr,10))/dfloat(integr/10)
        else
          fac=1.d0/dfloat(integr)
        endif
        do ii=1,3
          v(ii)=v(ii)*fac
        enddo
        i=i-1
        goto 90
      endif
      if(t(i).ne.'('.and.t(i).ne.'[') then
        if(maswrk) write(iun,*)
     .  '*** vec must be "X","Y","Z","D" or start with "(" or "["'
        call ending
      endif
      if(t(i).eq.'(') then
        ibracket=1
      else if (t(i).eq.'[') then
        ibracket=2
      endif
      do m=1,3
        i=i+1
        call parssr(t,nt,i,r)
        if(m.lt.3.and.t(i).ne.',') then
          if(maswrk) write(iun,*)
     .     '*** vec components must be seperated by ","'
          call ending
        endif
        if (ibracket.eq.1) then
          do jj=1,3
            v(jj)=v(jj)+rbas(jj,m)*r
          enddo
        else if (ibracket.eq.2) then
          v(m)=r
        endif
      enddo
      if(ibracket.eq.1.and.t(i).ne.')') then
        if(maswrk) write(iun,*)'*** ")" expected to close vec'
        call ending
      endif
      if(ibracket.eq.2.and.t(i).ne.']') then
        if(maswrk) write(iun,*)'*** "]" expected to close vec'
        call ending
      endif
  90  i=i+1
      end
