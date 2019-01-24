      subroutine sgroup(gen,gen_j,agen,ngen,g,g_j,ag0,ng,ngmx,lmax,jmax,
     &	                h,h_j,hh,hh_j,sig,sig_j,su1,su1_j,ifi)
c  Sets up space group given generators (gen,ag0en).
c  Operations are defined as  (g,a)(v):=g*v+a
c  where g is (3x3) matrix, a=vector.
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer :: ifi(*),ifih,ifihh,ifisig,ngmx,lmax,jmax,l,m,m1,iwig,
     &	         ist,ii,iwig0,jj,ng,igen,ngen,ig,itry,iord,nnow,j,ip0,i,
     &           k,n2,m2,isand,nnew,n
      logical spgeql
      real*8 :: gen(maxwig,*),g(maxwig,ngmx),h(maxwig),hh(maxwig),
     &          sig(maxwig),asig(3),agen(3,*),ag0(3,ngmx),
     &          ah(3),ahh(3),su1(-lmax:lmax,-lmax:lmax),ss
	complex*16 gen_j(maxwig_j,*),g_j(maxwig_j,ngmx),h_j(maxwig_j),
     &	         hh_j(maxwig_j),sig_j(maxwig_j),su1_j(jmax+1,jmax+1)
      do l=0,lmax
        do m=-l,l
          do m1=-l,l
            iwig=l*(2*l-1)*(2*l+1)/3+(2*l+1)*(l+m1)+l+m+1
            if(m.eq.m1)then
               g(iwig,1)=1.d0
            else
               g(iwig,1)=0.d0
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
                g_j(iwig,1)=dcmplx(ss,0.d0)
              enddo
	      enddo
          enddo
        enddo
	endif
      call dcopy(3,0.d0,0,ag0(1,1),1)
      ifia(1)=1
      ng=1
      do igen=1,ngen
        call spgcop(gen(1,igen),gen_j(1,igen),agen(1,igen),sig,sig_j,
     &              asig,ifi(igen),ifisig)
c ----- extend the group by all products with sig ----
        do ig=1,ng
          if(spgeql(g(1,ig),ag0(1,ig),sig,asig,ifia(ig),ifisig,gbas,
     &              maxwig)) then
            if(maswrk) write(iun,650) igen,ig
  650      format(' generator',i3,'  is already in group as element',i3)
            goto 80
          endif
        enddo
c  determine order
        call spgcop(sig,sig_j,asig,h,h_j,ah,ifisig,ifih)
        do itry=1,100
          iord=itry
          if(spgeql(h,ah,g(1,1),ag0(1,1),ifih,1,gbas,maxwig)) goto 2
          call spgprd(sig,sig_j,asig,h,h_j,ah,ifisig,ifih,hh,hh_j,ahh,
     &		          ifihh,lmax,jmax,su1,su1_j)
          call spgcop(hh,hh_j,ahh,h,h_j,ah,ifihh,ifih)
        enddo
c  products of type:  g1 sig**p g2
2       nnow=ng
        if(maswrk) write(iun,677) igen,iord
677     format(' Order of gen',i3,'  is',i5)
        n2=0
        do j=1,ng
          call spgcop(g(1,j),g_j(1,j),ag0(1,j),h,h_j,ah,ifia(j),ifih)
          do ip0=1,iord-1
            call spgprd(sig,sig_j,asig,h,h_j,ah,ifisig,ifih,hh,hh_j,
     & 		            ahh,ifihh,lmax,jmax,su1,su1_j)
            call spgcop(hh,hh_j,ahh,h,h_j,ah,ifihh,ifih)
            do i=1,ng
              call spgprd(g(1,i),g_j(1,i),ag0(1,i),h,h_j,ah,ifia(i),
     &	         		  ifih,hh,hh_j,ahh,ifihh,lmax,jmax,su1,su1_j)
              do k=1,nnow
                if( spgeql(g(1,k),ag0(1,k),hh,ahh,ifia(k),ifihh,
     &                          gbas,maxwig)) goto 11
              enddo
  13          nnow=nnow+1
              if(nnow.gt.ngmx) goto 99
              call spgcop(hh,hh_j,ahh,g(1,nnow),g_j(1,nnow),
     &			          ag0(1,nnow),ifihh,ifia(nnow))
11            continue
            enddo
          enddo
          if(j.eq.1) n2=nnow
        enddo
c   products with more than one sandwiched sigma-factor
        m1=ng+1
        m2=nnow
        do isand=2,50
          nnew=0
          do n=ng+1,n2
            do m=m1,m2
              call spgprd(g(1,n),g_j(1,n),ag0(1,n),g(1,m),g_j(1,m),
     &			           ag0(1,m),ifia(n),ifia(m),h,h_j,ah,ifih,lmax,
     &                     jmax,su1,su1_j)
              do k=1,nnow
                if(spgeql(g(1,k),ag0(1,k),h,ah,ifia(k),ifih,gbas,
     &                      maxwig)) goto 21
              enddo
              nnew=nnew+1
              nnow=nnow+1
              if(nnow.gt.ngmx) goto 99
              call spgcop(h,h_j,ah,g(1,nnow),g_j(1,nnow),ag0(1,nnow),
     &			          ifih,ifia(nnow))
21            continue
            enddo
          enddo
          m1=m2+1
          m2=nnow
          if(nnew.eq.0) goto 25
        enddo
25      continue
        ng=nnow
80      continue
      enddo
      return
  99  if(maswrk) write(iun,*)
     &   'SGROUP : **** ng greater than     ngmx=',ngmx
      call ending
      end
