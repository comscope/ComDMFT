      subroutine get_str_con(sc,pnt0)

c     make matrix of all structure constants for one k-point, e = 0
!     nonlocal part of MT-MT coulmb matrix

      use atom_mod
      use solid_mod
      use units_mod

      implicit none
      complex*16, intent(out) :: sc(nlb2,nlb2)
      integer :: nlmi0,lmaxst,nlmst,id,jb,lmaxb,ilb0,ib,
     &  ip0,ipa,m,lmaxa,ila0
      double precision :: sig,scalp,dlatm,pnt0(3)
      complex*16 eiphi
      complex*16, allocatable :: dl(:),b1(:,:),b2(:,:)

      nlmi0=limlpb
      lmaxst=maxpb+maxpb
      nlmst=(lmaxst+1)**2
      allocate(dl(nlmst))
      allocate(b1(nlmi0,nlmi0))
      allocate(b2(nlmi0,nlmi0))

c     ------ start loop over connecting vectors -----
      do id=1,ndtau
c     ---  Get H(l") for pair id with dtau=r(ib)-r(jb) -------------------
        call rcnke0(pnt0,dtau(1,id),lmaxst,dl)
c     ------ Put H(l") into SS(l,l') -------------------------------------
        call cansl0(b1,nlmi0,nlmi0,dl)
c     ----------- SS for pair with -dtau --------------------------------
        call xxmste0(b1,b2,nlmi0)
c     --- sort into arrays b,bb. jb=center, ib=where expanded ----
        do jb=1,natom
          lmaxb=lmpb(is(jb))
          ilb0=iob2(jb)-1
          do ib=1,natom
            ip0=ipdtau(ib,jb)   !! if>0 means ib <- jb, otherwise ib -> jb
            ipa=iabs(ip0)
            if(ipa.eq.id) then
              sig=1.d0
              if(ip0.lt.0) sig=-1.d0
              scalp=0.d0
              do m=1,3
                dlatm=tau(m,ib)-tau(m,jb)-sig*dtau(m,ipa) ! lat.vec.ne.0
                scalp=scalp+2.d0*pi*dlatm*pnt0(m)
              enddo
              eiphi=dcmplx(dcos(scalp),dsin(scalp)) !! phase factor
              lmaxa=lmpb(is(ib))
              ila0=iob2(ib)-1
c     -------------- to just multiply on possible phase factor -------------
              if(ip0.gt.0) then
                call yymste(b1,lmaxa,lmaxb,nlmi0,nlmi0,sc,nlb2,nlb2,
     &            ila0,ilb0,eiphi)
              else
                call yymste(b2,lmaxa,lmaxb,nlmi0,nlmi0,sc,nlb2,nlb2,
     &            ila0,ilb0,eiphi)
              endif
            endif
          enddo                 !! over jb
        enddo                   !! over ib
      enddo                     !! over id
      deallocate(dl,b1,b2)
      end
