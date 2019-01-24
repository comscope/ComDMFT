

      subroutine integral_band_pw_band_wan(ff,znew,anew,znew1,anew1,
     $  b_pw_b,k,pn,kq,kq_shift)

      use comwann_mod

      implicit none
      integer, intent(in) :: k,kq,kq_shift(3)
      double precision, intent(in) :: pn(3)
      complex*16, intent(in) ::
     $  znew(nfun,num_bands),
     &  anew(nbasmpw,num_bands),
     &  znew1(nfun,num_bands),
     &  anew1(nbasmpw,num_bands),
     &  ff(maxel,maxel,natom)
      complex*16, intent(out) :: b_pw_b(num_bands,num_bands)
      integer :: j,iatom,isort,k0,kq0,ndimb,ind0,
     &  jbas,jnd,ia,ib,ic,ibas,ind,ia1,
     &  ib1,ic1,ibas1,ibas0,irl
!      double precision :: v(3),det,det3
      complex*16, allocatable :: s(:,:),tmp(:,:),tmp1(:,:)

      allocate(s(nbndtm,nbndtm))
      allocate(tmp(nbndtm,num_bands))
      allocate(tmp1(maxel,num_bands))
      s=0.0d0

      b_pw_b=(0.d0,0.d0)

      k0=i_kref(k)
      kq0=i_kref(kq)
c     ----------------- MT contribution ----------------------------------
      do iatom=1,natom
        isort=is(iatom)
        ndimb=lfunm(isort)
        ind0=io_lem(iatom)
        tmp1=0.0d0
        call zgemm('n','n',ndimb,num_bands,ndimb,
     &    (1.d0,0.d0),ff(1,1,iatom),maxel,
     &    znew1(ind0,1),nfun,(0.d0,0.d0),tmp1,maxel)
        call zgemm('c','n',num_bands,num_bands,ndimb,
     &    (1.d0,0.d0),znew(ind0,1),nfun,
     &    tmp1,maxel,(1.d0,0.d0),b_pw_b,
     &    num_bands)
      enddo                     !! over iatom
c     ----------------- Interstitial contribution ------------------------
      do jbas=1,nbask(k0)/nrel
        jnd=indgb(jbas,k)
        ia=kq_shift(1)+igbs(1,jnd)
        ib=kq_shift(2)+igbs(2,jnd)
        ic=kq_shift(3)+igbs(3,jnd)
        do ibas=1,nbask(kq0)/nrel
          ind=indgb(ibas,kq)
          ia1=-ia+igbs(1,ind)
          ib1=-ib+igbs(2,ind)
          ic1=-ic+igbs(3,ind)
          j=indplw(ia1,ib1,ic1)
          if(complex_ro) then
            s(jbas,ibas)=dcmplx(sovr(j),sovi(j))
          else
            s(jbas,ibas)=dcmplx(sovr(j),0.d0)
          endif
        enddo
      enddo
      do irl=1,nrel
        ibas1=nbask(kq0)/nrel*(irl-1)
        ibas0=nbask(k0)/nrel*(irl-1)
        tmp=0.0d0
        call zgemm('n','n',nbask(k0)/nrel,num_bands,nbask(kq0)/nrel,
     &    (1.d0,0.d0),s,nbndtm,anew1(ibas1+1,1),nbasmpw,
     &    (0.d0,0.d0),tmp,nbndtm)
        call zgemm('c','n',num_bands,num_bands,
     &    nbask(k0)/nrel,(1.d0,0.d0),anew(ibas0+1,1),nbasmpw,
     &    tmp,nbndtm,(1.d0,0.d0),b_pw_b,
     &    num_bands)
      enddo
      deallocate(s,tmp,tmp1)
      end

