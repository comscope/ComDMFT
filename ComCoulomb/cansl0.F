      subroutine cansl0(s,nlma,nlma0,hl)

c     makes lattice structure constants. make hl
c     in rcnke0 but note that there dr is in alat units, q in tpiba units.
c     exp. theorem is:  k(k,r-dr) = sum(m) s(m,k,dr)*j(m,r)
c     functions are real hankels and bessels without extra factors. E=0

      use manager_mod
      use parallel_mod
      implicit none
      integer :: nlma,nlma0,lmax,l,mlm,lm,lget,klm,lk,ii,icg,llm,ll0
      double precision :: sig(0:39),fpi,fac2
      complex*16 :: s(nlma0,nlma0),hl(*),sum
      fpi=16.d0*datan(1.d0)
      lmax=lget(nlma)+lget(nlma)
      if(lmax.gt.39) then
        if(maswrk) write(iun,*)'*** change dimensions in cansle0'
c        write(179,*) "calling ending in cansl0.f"
        call ending
      endif
      sig(0)=1.d0
      do 1 l=1,lmax
 1      sig(l)=-sig(l-1)
c     ---------------------------------------
      do mlm=1,nlma
        lm=lget(mlm)
        do klm=1,nlma
          lk=lget(klm)
          sum=(0.d0,0.d0)
          ii=max0(mlm,klm)
          ii=(ii*(ii-1))/2+min0(mlm,klm)
          do icg=indxcg(ii),indxcg(ii+1)-1
            llm=jcg(icg)
            ll0=lget(llm)
            if(ll0.eq.lm+lk) sum=sum+cg(icg)*hl(llm)
          enddo
          s(mlm,klm)=2.d0*fpi*sum*sig(lk)/fac2(lm+lm+1)/fac2(lk+lk+1)
        enddo
      enddo
      end
      
