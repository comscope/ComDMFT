
      subroutine fifi_j_prepare_onepoint(fifi_j,pn,ispin)

      use comwann_mod
      implicit none
      integer,intent(in) :: ispin
      double precision, intent(in) :: pn(3)
      double precision, intent(out) ::
     $  fifi_j(maxlfun,maxlfun,0:2*maxb,
     &  nsort)
      integer :: isort,ir,l,mt,mt1,lf,lf1
      double precision :: pi4,qq,rv,dqdall,q(3),c2,pp,pi2a,
     $  bess(0:2*maxb,0:maxnrad),bessd(0:2*maxb),
     $  work(0:maxnrad),work1(0:maxnrad)

      fifi_j=0.0d0

      pi4=4.d0*pi
      pi2a=(pi+pi)/par
!      sqpi4=sqrt(pi4)
      c2=clight**2

      q=pn
      qq=sqrt(dot_product(q,q))*pi2a ! bvec
      
      do isort=1,nsort
        do ir=0,nrad(isort)
          rv=r(ir,isort)
          call BESSR(qq,rv,bess(0,ir),bessd,2*lmb(isort))
          rv=rv*rv
          work(ir)=rv*dr(ir,isort)
        enddo
        do l=0,2*lmb(isort)
          do lf1=1,lfun(isort)
            mt1=ind_wf(lf1,isort)
            do lf=1,lfun(isort)
              mt=ind_wf(lf,isort)
              do ir=0,nrad(isort)
                pp=gfun(mt+ir,ispin)*gfun(mt1+ir,ispin)
                if(irel>=1) pp=pp+gfund(mt+ir,ispin)
     &            *gfund(mt1+ir,ispin)/c2
                work1(ir)=work(ir)*bess(l,ir)*pp
              enddo
              fifi_j(lf,lf1,l,isort)=
     &          pi4*dqdall(h(isort),work1,nrad(isort))
            enddo               !! over lf
          enddo                 !! over lf1
        enddo                   !! over l
      enddo                     !! over iatom
      end      

