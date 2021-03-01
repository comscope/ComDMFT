      subroutine symg_loc(gloc)
      use atom_mod
      use manager_mod
      use solid_mod
      implicit none
      integer :: iatom,isort,ig,jatom,j,l,lp1,ie,m,lm,km,lm0,li,
     &  ist,ii,jj,mj,i,in,jn
      complex*16 :: cf,gloc(maxel,maxel,natom)
      complex*16, allocatable :: y0j(:),t0j(:),g_tmp(:,:,:),pa1(:,:),
     &  y0(:),t0(:)
      allocate(y0j(2*maxb+2))
      allocate(t0j(2*maxb+2))
      allocate(g_tmp(maxel,maxel,natom))
      allocate(pa1(maxel,maxel))
      allocate(y0(2*maxb+1))
      allocate(t0(2*maxb+1))
      cf=dcmplx(1.d0/dfloat(ngroup),0.d0)
      g_tmp=gloc
      gloc=(0.d0,0.d0)
      do iatom=1,natom
        isort=is(iatom)
        do ig=1,ngroup
          jatom=ip(iatom,ig)
          do j=1,lfunm(isort)
            if(irel.le.1) then
              do l=0,lmb(isort)
                lp1=l+1
                do ie=1,ntle(l,isort)
                  in=1
                  if(augm(ie,l,isort)/='LOC') in=2
                  do jn=1,in
                    do m=-l,l
                      lm=l*(l+1)+m+1
                      km=indbasa(jn,ie,lm,isort)
                      lm0=l+m+1
                      t0(lm0)=g_tmp(km,j,jatom)
                    enddo
                    call rotate1c(t0,y0,l,u(1,ig),1)
                    do m=-l,l
                      lm=l*(l+1)+m+1
                      km=indbasa(jn,ie,lm,isort)
                      lm0=l+m+1
                      pa1(km,j)=y0(lm0)
                    enddo
                  enddo
                enddo
              enddo             !! over l
            else if(irel.eq.2) then
              li=0
              do l=0,lmb(isort)
                ist=-1
                if(l.eq.0) ist=1
                do ii=ist,1,2
                  li=li+1
                  jj=l+l+ii
                  do ie=1,ntle(l,isort)
                    in=1
                    if(augm(ie,l,isort)/='LOC') in=2
                    do jn=1,in
                      lm0=0
                      do mj=-jj,jj,2
                        lm0=lm0+1
                        call getlimj(lm,l,ii,mj,li,1)
                        km=indbasa(jn,ie,lm,isort)
                        t0j(lm0)=g_tmp(km,j,jatom)
                      enddo
                      call rotate2(t0j,y0j,l,ii,uj(1,ig),1,1)
                      do mj=-jj,jj,2
                        call getlimj(lm,l,ii,mj,li,1)
                        km=indbasa(jn,ie,lm,isort)
                        lm0=(jj+mj)/2+1
                        pa1(km,j)=y0j(lm0)
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            endif
          enddo                 !! over j
          do i=1,lfunm(isort)
            if(irel.le.1) then
              do l=0,lmb(isort)
                lp1=l+1
                do ie=1,ntle(l,isort)
                  in=1
                  if(augm(ie,l,isort)/='LOC') in=2
                  do jn=1,in
                    do m=-l,l
                      lm=l*(l+1)+m+1
                      km=indbasa(jn,ie,lm,isort)
                      lm0=l+m+1
                      t0(lm0)=pa1(i,km)
                    enddo
                    call rotate1c(t0,y0,l,u(1,ig),1)
                    do m=-l,l
                      lm=l*(l+1)+m+1
                      km=indbasa(jn,ie,lm,isort)
                      lm0=l+m+1
                      gloc(i,km,iatom)=gloc(i,km,iatom)+y0(lm0)
                    enddo
                  enddo
                enddo
              enddo
            else if(irel.eq.2) then
              li=0
              do l=0,lmb(isort)
                ist=-1
                if(l.eq.0) ist=1
                do ii=ist,1,2
                  li=li+1
                  jj=l+l+ii
                  do ie=1,ntle(l,isort)
                    in=1
                    if(augm(ie,l,isort)/='LOC') in=2
                    do jn=1,in
                      lm0=0
                      do mj=-jj,jj,2
                        lm0=lm0+1
                        call getlimj(lm,l,ii,mj,li,1)
                        km=indbasa(jn,ie,lm,isort)
                        t0j(lm0)=pa1(i,km)
                      enddo
                      call rotate2(t0j,y0j,l,ii,uj(1,ig),1,2)
                      do mj=-jj,jj,2
                        call getlimj(lm,l,ii,mj,li,1)
                        km=indbasa(jn,ie,lm,isort)
                        lm0=(jj+mj)/2+1
                        gloc(i,km,iatom)=gloc(i,km,iatom)+y0j(lm0)
                      enddo
                    enddo
                  enddo         !! over ie
                enddo           !! over ii
              enddo             !! over l
            endif
          enddo                 !! over i
        enddo                   !! over ig
      enddo                     !! over iatom
      gloc=cf*gloc
      deallocate(y0j,t0j,g_tmp,pa1,y0,t0)
      end
