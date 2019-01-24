      subroutine integral_band_pw_band_mt_wan(fifi_j,ff,pn)
      use comwann_mod
      implicit none
      include 'mpif.h'
      double precision, intent(in) ::
     $  fifi_j(maxlfun,maxlfun,0:2*maxb,nsort),pn(3)
      complex*16, intent(out) :: ff(maxel,maxel,natom)

      integer :: iatom,isort,ndimb,lm1,lm2,l1,l2,ii,
     &  iii,indx,icg1,icg2,i2,mj1,mj2,i1,li1,li2,icg,lmi,lget,
     &  l_i,nm1,nm2,ie1,km1,lf1,ie2,km2,lf2,in1,
     &  jn2,in2,jn1
      double precision :: fas,q(3),qq,sqpi4,pi2
      complex*16 :: c1,c2
      double precision, allocatable :: ylm(:)

!      write(180,*)'indxcgr',indxcgr
!      call mpi_barrier(mpi_comm_world,mpierr)       
 
      allocate(ylm((2*maxb+1)**2))
      pi2=pi+pi
      sqpi4=sqrt(4.d0*pi)
      q=pn
      qq=sqrt(dot_product(q,q))

      if(qq.gt.1.d-9) then
        call sphharm(q(1)/qq,q(2)/qq,q(3)/qq,2*maxb,ylm)
      else
        ylm=0.d0
        ylm(1)=1.d0/sqpi4
      endif
      ff=(0.d0,0.d0)
      
      do iatom=1,natom
        isort=is(iatom)
!        nl=lfunm(isort)
        ndimb=nrel*(lmb(isort)+1)**2
        fas=pi2*dot_product(pn,tau(:,iatom))
        c1=exp(dcmplx(0.d0,-fas))
        do lm1=1,ndimb
          do lm2=1,ndimb
            if(irel<=1) then
              l1=lget(lm1)
              l2=lget(lm2)
              ii = max0(lm1,lm2)
              iii = min0(lm1,lm2)
              indx = (ii*(ii-1))/2 + iii
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1) - 1 
            else
              call getlimj(lm2,l2,i2,mj2,li2,0)
              call getlimj(lm1,l1,i1,mj1,li1,0)
!              write(179,*) 'lm1',lm1,'lm2',lm2
              icg1=indxcgr(lm1,lm2,1,0)
              icg2=indxcgr(lm1,lm2,2,0)
            endif
            do icg = icg1, icg2
              if(irel<=1) then
                lmi = jcg(icg)
                l_i=lget(lmi)
                c2=c1*dcmplx(0.d0,-1.d0)**l_i*ylm(lmi)*cg(icg)
              else
                lmi=jcgr(icg)
                l_i=lget(lmi)
                c2=c1*dcmplx(0.d0,-1.d0)**l_i*ylm(lmi)*cgr(icg)
              endif
              nm1=ntle(l1,isort)
              nm2=ntle(l2,isort)
              do ie1=1,nm1
                in1=1
                if(augm(ie1,l1,isort)/='LOC') in1=2
                do jn1=1,in1
                  km1=indbasa(jn1,ie1,lm1,isort)
                  lf1=lf_isz(km1,isort)
                  do ie2=1,nm2
                    in2=1
                    if(augm(ie2,l2,isort)/='LOC') in2=2
                    do jn2=1,in2
                      km2=indbasa(jn2,ie2,lm2,isort)
                      lf2=lf_isz(km2,isort)
                      ff(km1,km2,iatom)=
     &                  ff(km1,km2,iatom)+c2
     &                  *fifi_j(lf1,lf2,l_i,isort)
                    enddo
                  enddo         ! over ie2
                enddo
              enddo
            enddo               ! over icg
          enddo                 ! over lm2
        enddo                   ! over lm1
      enddo                     !!! over iatom
      deallocate(ylm)
      end

