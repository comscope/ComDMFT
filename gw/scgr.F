      subroutine scgr(icg,key)
c     key=0 - first pass
c     key=1 - second pass
	use atom_mod
	use manager_mod
      implicit none
	integer, intent(in) :: key
	integer, intent(out) :: icg
	integer :: mcase,lm,l,ist,i,jj,mj,lm2,li,l1m1,l1,is1,i1,jj1,mj1,
     &	         icase,lm1,l0,lget,m0
	real*8 :: angr,angi
      mcase=0
      if(magn.eq.2) mcase=3
      lm=0
      icg=0
      do l=0,maxb+1
        ist=-1
        if(l.eq.0) ist=1
        do i=ist,1,2
          jj=l+l+i
          do mj=-jj,jj,2
            lm=lm+1
            call getlimj(lm2,l+i,-i,mj,li,1)
            lm_tild(lm)=lm2
            l1m1=0
            do l1=0,maxb+1
              is1=-1
              if(l1.eq.0) is1=1
              do i1=is1,1,2
                jj1=l1+l1+i1
                do mj1=-jj1,jj1,2
                  l1m1=l1m1+1
                  do icase=0,mcase
                    indxcgr(lm,l1m1,1,icase)=icg+1
                    do lm1=1,limlpb
                      l0=lget(lm1)
                      m0=lm1-l0*(l0+1)-1
                      call angle(l0,m0,l,i,mj,l1,i1,mj1,icase,angr,angi)
                      if(abs(angr).gt.1.d-12.or.abs(angi).gt.1.d-12)then
                        icg=icg+1
                        jcgr(icg)=lm1
                        if(key.eq.1) cgr(icg)=dcmplx(angr,angi)
                      endif
                    enddo
                    indxcgr(lm,l1m1,2,icase)=icg
                  enddo  !!! over icase
                enddo  !!! over mj1
              enddo  !!! over i1
            enddo  !!! over l1
          enddo  !!! over mj
        enddo  !!! over i
      enddo  !!! over l
      end