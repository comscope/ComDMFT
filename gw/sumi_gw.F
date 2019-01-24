      subroutine sumi_gw(correlation,ind_k,k,ispin)
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use lda_dmft_mod
      implicit none
      logical, intent(in) :: correlation
      integer :: ibas,jbas,ind,jnd,ia,ib,ic,j,k,ispin,igb,ii,jgb,ipp,ih,
     &  jp,jh,i,ndim,ind_k,nbas,ibas0,irl,n0,jbas0
      real*8 :: cr(3),gki(3),gkj(3),wei,pi2,pi2a,c05,c2,g2i,gi,eni,
     &  c2i,const1,anormi,g2j,gj,enj,c2j,anormj,const,const2,
     &  scal,ddot,ann,sr,si,womr(3),womi(3),relr,reli
      complex*16 :: cc
      complex*16, allocatable :: tmp(:,:),aga(:,:)
      n0=nbndf
      ndim=n_bnd(k,ispin)
      wei=wgt(k)*2.d0/dfloat(nrel*nspin)
      if(irel<=1) then
        nbas=nbask(k)
        allocate(tmp(nbas,ndim),aga(nbas,nbas))
        if(.not.correlation) then
          if (is_lda_dmft) then
            call zgemm('n','n',nbas,ndim,ndim,(1.d0,0.d0),
     &        ev_bnd(1,1,ind_k,ispin),nbasmpw,
     &        g_full_0(1,1,ind_k,ispin),n0,(0.d0,0.d0),tmp,nbas)            
          else
            do j=1,ndim
              do i=1,nbas
                tmp(i,j)=ev_bnd(i,j,ind_k,ispin)
     $            *g_full_00(j,ind_k,ispin)
              enddo
            enddo
          endif
        else
          call zgemm('n','n',nbas,ndim,ndim,(1.d0,0.d0),
     &      ev_bnd(1,1,ind_k,ispin),nbasmpw,
     &      g_full_0(1,1,ind_k,ispin),n0,(0.d0,0.d0),tmp,nbas)
        endif
        call zgemm('n','c',nbas,nbas,ndim,(1.d0,0.d0),tmp,nbas,
     &    ev_bnd(1,1,ind_k,ispin),nbasmpw,(0.d0,0.d0),aga,
     &    nbas)
        deallocate(tmp)
        aga=aga/amega
c     --------------------------------------------------------------------
        do ibas=1,nbas
          ind=indgb(ibas,k)
          do jbas=1,nbas
            jnd=indgb(jbas,k)
            ia=igbs(1,ind)-igbs(1,jnd)
            ib=igbs(2,ind)-igbs(2,jnd)
            ic=igbs(3,ind)-igbs(3,jnd)
            j=indplw(ia,ib,ic)
            if(j>nplwro) cycle
            rointr(j,ispin)=rointr(j,ispin)+wei*aga(ibas,jbas)
          enddo                 !!! over jbas
        enddo                   !!! over ibas
        deallocate(aga)
      else if(irel==2) then
        if(.not.rel_interst) then
          nbas=nbask(k)/2
          allocate(tmp(nbas,ndim),aga(nbas,nbas))
          do irl=1,2
            ibas0=0
            if(irl==2) ibas0=nbas
            if(.not.correlation) then
              if (is_lda_dmft) then
                call zgemm('n','n',nbas,ndim,ndim,(1.d0,0.d0),
     &            ev_bnd(ibas0+1,1,ind_k,1),nbasmpw,
     &            g_full_0(1,1,ind_k,1),n0,(0.d0,0.d0),tmp,nbas)                
              else
                do j=1,ndim
                  do i=1,nbas
                    tmp(i,j)=ev_bnd(ibas0+i,j,ind_k,ispin)
     &                *g_full_00(j,ind_k,ispin)
                  enddo
                enddo
              endif
            else
              call zgemm('n','n',nbas,ndim,ndim,(1.d0,0.d0),
     &          ev_bnd(ibas0+1,1,ind_k,1),nbasmpw,
     &          g_full_0(1,1,ind_k,1),n0,(0.d0,0.d0),tmp,nbas)
            endif
            call zgemm('n','c',nbas,nbas,ndim,(1.d0,0.d0),tmp,nbas,
     &        ev_bnd(ibas0+1,1,ind_k,1),nbasmpw,(0.d0,0.d0),
     &        aga,nbas)
            aga=aga/amega
c     --------------------------------------------------------------------
            do ibas=1,nbas
              ind=indgb(ibas,k)
              do jbas=1,nbas
                jnd=indgb(jbas,k)
                ia=igbs(1,ind)-igbs(1,jnd)
                ib=igbs(2,ind)-igbs(2,jnd)
                ic=igbs(3,ind)-igbs(3,jnd)
                j=indplw(ia,ib,ic)
                if(j>nplwro) cycle
                cc=wei*aga(ibas,jbas)
                rointr(j,1)=rointr(j,1)+cc
                if(magn.eq.2) then
                  spintr(3,j)=spintr(3,j)+(-1)**irl*cc
                endif
              enddo             !!! over jbas
            enddo               !!! over ibas
          enddo                 !! over irl
          if(magn.eq.2) then
            ibas0=0
            jbas0=nbas
            if(.not.correlation) then
              if (is_lda_dmft) then
                call zgemm('n','n',nbas,ndim,ndim,(1.d0,0.d0),
     &            ev_bnd(ibas0+1,1,ind_k,1),nbasmpw,
     &            g_full_0(1,1,ind_k,1),n0,(0.d0,0.d0),tmp,nbas)                
              else
                do j=1,ndim
                  do i=1,nbas
                    tmp(i,j)=ev_bnd(ibas0+i,j,ind_k,ispin)
     &                *g_full_00(j,ind_k,ispin)
                  enddo
                enddo
              endif
            else
              call zgemm('n','n',nbas,ndim,ndim,(1.d0,0.d0),
     &          ev_bnd(ibas0+1,1,ind_k,1),nbasmpw,
     &          g_full_0(1,1,ind_k,1),n0,(0.d0,0.d0),tmp,nbas)
            endif
            call zgemm('n','c',nbas,nbas,ndim,(1.d0,0.d0),tmp,nbas,
     &        ev_bnd(jbas0+1,1,ind_k,1),nbasmpw,(0.d0,0.d0),
     &        aga,nbas)
            aga=aga/amega
c     --------------------------------------------------------------------
            do ibas=1,nbas
              ind=indgb(ibas,k)
              do jbas=1,nbas
                jnd=indgb(jbas,k)
                ia=igbs(1,ind)-igbs(1,jnd)
                ib=igbs(2,ind)-igbs(2,jnd)
                ic=igbs(3,ind)-igbs(3,jnd)
                j=indplw(ia,ib,ic)
                if(j>nplwro) cycle
                cc=wei*aga(ibas,jbas)
                spintr(1,j)=spintr(1,j)+cc
                spintr(2,j)=spintr(1,j)+(0.d0,1.d0)*cc
                cc=wei*conjg(aga(jbas,ibas))
                spintr(1,j)=spintr(1,j)+cc
                spintr(2,j)=spintr(1,j)-(0.d0,1.d0)*cc
              enddo             !!! over jbas
            enddo               !!! over ibas
          endif
          deallocate(aga,tmp)
        else if(rel_interst) then
          nbas=nbask(k)
          allocate(tmp(nbas,ndim),aga(nbas,nbas))
          if(.not.correlation) then
            if (is_lda_dmft) then
              call zgemm('n','n',nbas,ndim,ndim,(1.d0,0.d0),
     &          ev_bnd(1,1,ind_k,1),nbasmpw,
     &          g_full_0(1,1,ind_k,1),n0,(0.d0,0.d0),tmp,nbas)
            else
              do j=1,ndim
                do i=1,nbas
                  tmp(i,j)=ev_bnd(i,j,ind_k,ispin)
     &              *g_full_00(j,ind_k,ispin)
                enddo
              enddo
            endif
          else
            call zgemm('n','n',nbas,ndim,ndim,(1.d0,0.d0),
     &        ev_bnd(1,1,ind_k,1),nbasmpw,
     &        g_full_0(1,1,ind_k,1),n0,(0.d0,0.d0),tmp,nbas)
          endif
          call zgemm('n','c',nbas,nbas,ndim,(1.d0,0.d0),tmp,nbas,
     &      ev_bnd(1,1,ind_k,ispin),nbasmpw,(0.d0,0.d0),aga,
     &      nbas)
          deallocate(tmp)
          aga=aga/amega
          pi2=2.d0*pi
          pi2a=pi2/par
          c05=0.5d0*clight
          c2=clight*clight
          do igb=1,nbask(k)/2
            ind=indgb(igb,k)
            g2i=0.d0
            do ii=1,3
              gki(ii)=(pnt(ii,k)+gbs(ii,ind))*pi2a
              g2i=g2i+gki(ii)**2
            enddo
            gi=sqrt(g2i)
            eni=c05*(sqrt(c2+(2.d0*gi)**2)-clight)
            c2i=c2+eni
            const1=clight/c2i
            anormi=c2i/sqrt(c2i*c2i+(clight*gi)**2)
            do jgb=1,nbask(k)/2
              jnd=indgb(jgb,k)
              ia=igbs(1,ind)-igbs(1,jnd)
              ib=igbs(2,ind)-igbs(2,jnd)
              ic=igbs(3,ind)-igbs(3,jnd)
              j=indplw(ia,ib,ic)
              if(j>nplwro) cycle
              g2j=0.d0
              do ii=1,3
                gkj(ii)=(pnt(ii,k)+gbs(ii,jnd))*pi2a
                g2j=g2j+gkj(ii)**2
              enddo
              gj=sqrt(g2j)
              enj=c05*(sqrt(c2+(2.d0*gj)**2)-clight)
              c2j=c2+enj
              anormj=c2j/dsqrt(c2j*c2j+(clight*gj)**2)
              const=c2/c2i/c2j
              const2=clight/c2j
              scal=ddot(3,gki,1,gkj,1)
              call cross(cr,gki,gkj)
              ann=anormi*anormj
              do ipp=1,2
                ih=(ipp-1)*nbask(k)/2+igb
                do jp=1,2
                  jh=(jp-1)*nbask(k)/2+jgb
                  sr=dreal(aga(ih,jh))
                  si=dimag(aga(ih,jh))
                  call setrel(cr,scal,const,ipp,jp,ann,gki,gkj,womr,
     &              womi,relr,reli,magn)
                  rointr(j,1)=rointr(j,1)+wei*(sr*relr-si*reli)
                  if(complex_ro) rointr(j,1)=rointr(j,1)
     &              +(0.d0,1.d0)*wei
     &              *(sr*reli+si*relr)
                  if(magn.eq.2) then
                    do i=1,3
                      spintr(i,j)=spintr(i,j)+wei
     &                  *dcmplx(sr*womr(i)-si*womi(i),
     &                  sr*womi(i)+si*womr(i))
                    enddo
                  endif
                enddo           !!! jp
              enddo             !!! ipp
            enddo               !!! over jgb
          enddo                 !!! over igb
          deallocate(aga)
        endif
      endif                     !!! over irel
      end
