      subroutine fourier(rhors,rhois,indstr,roin,key,igbsf0,
     &                   gbs0,ia,ib,jreal,jimag)
c     key=0 - from direct to reciprocal
c     key=1 - from reciprocal to direct
	use solid_mod
      IMPLICIT none
      INTEGER indstr(0:nstar),igbsf0(nplwro),key,ia,ib,isi,istar,j,ifft
      logical jreal,jimag
      real*8 :: rhors(mmm),rhois(mmm),gbs0(3,nplwro),rh,gbb,anorm
      complex*16 :: roin(*)
      isi=(-1)**(key+1)
      if(key.eq.1) then
        call dcopy(mmm,0.d0,0,rhors,1)
        call dcopy(mmm,0.d0,0,rhois,1)
        do istar=1,nstar
          do j=indstr(istar-1)+1,indstr(istar)
            ifft=igbsf0(j)
            if(jreal) rhors(ifft)=real(roin(j))
            if(jimag) rhois(ifft)=imag(roin(j))
            if(ia.ne.0.and.ib.eq.0) then
              rh=-rhois(ifft)*gbs0(ia,j)
              rhois(ifft)=rhors(ifft)*gbs0(ia,j)
              rhors(ifft)=rh
            endif
            if(ia.ne.0.and.ib.ne.0) then
              gbb=-gbs0(ia,j)*gbs0(ib,j)
              rhors(ifft)=rhors(ifft)*gbb
              rhois(ifft)=rhois(ifft)*gbb
            endif
          enddo
        enddo
      endif
      call cft(rhors,rhois,mmm,mdiv(1),mdiv(1),isi)
      call cft(rhors,rhois,mmm,mdiv(2),mdiv(1)*mdiv(2),isi)
      call cft(rhors,rhois,mmm,mdiv(3),mmm,isi)
C  **** SET PSEUDODENSITY WITHIN SPHERE GMAX ***
      if(key.eq.0) then
        ANORM=1.D0/DBLE(mmm)
        do j=1,nplwro
          ifft=igbsf0(j)
          if(jreal) roin(j)=roin(j)+anorm*rhors(ifft)
          if(jimag) roin(j)=roin(j)+(0.d0,1.d0)*anorm*rhois(ifft)
        enddo
      endif
      END