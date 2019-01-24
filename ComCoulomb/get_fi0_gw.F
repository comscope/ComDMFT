      subroutine get_fi0_gw

      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer :: ispin,isort,le1,mt1,mt,le,ir,i,mtpb,l,jspin
      double precision :: sqpi8,c2,dqdall
      double precision, allocatable :: pp(:),pp1(:),qq(:),qq1(:),work(:)
      allocate(pp(0:maxnrad),pp1(0:maxnrad),qq(0:maxnrad),
     $  qq1(0:maxnrad),work(0:maxnrad))


c$$$      matrix element calculation in the MT region

      c2=clight**2
      do jspin=1,nspin
        do ispin=1,nspin
          do isort=1,nsort
            do le1=1,lfun(isort)
              mt1=ind_wf(le1,isort)
              do ir=0,nrad(isort)
                pp1(ir)=gfun(mt1+ir,jspin)*r(ir,isort)
                if(irel>=1) qq1(ir)=gfund(mt1+ir,jspin)*r(ir,isort)
              enddo
              do le=1,lfun(isort)
                mt=ind_wf(le,isort)
                do ir=0,nrad(isort)
                  pp(ir)=gfun(mt+ir,ispin)*r(ir,isort)
                  if(irel>=1) qq(ir)=gfund(mt+ir,ispin)*r(ir,isort)
                enddo
                do i=1,lfun_pb(isort)
                  mtpb=ind_pb(i,isort)
                  do ir=0,nrad(isort)
                    work(ir)=pp(ir)*pp1(ir)
                    if(irel>=1) work(ir)=work(ir)+qq(ir)*qq1(ir)/c2
                    work(ir)=work(ir)*fun_pb(ir+mtpb)*dr(ir,isort)
                  enddo
                  fi0(le,le1,i,isort,ispin,jspin)=dqdall(h(isort),work,
     $              nrad(isort)) ! 
                enddo		!! over i
              enddo		!! over le
            enddo		!! over le1
          enddo                 !! over isort
        enddo			!! over ispin
      enddo			!! over jspin
c     -------- Setting Mi0 --- VII-13-5 -----------------------------------------

      sqpi8=sqrt(8.d0*pi)
      do isort=1,nsort
        do i=1,lfun_pb(isort)
          mtpb=ind_pb(i,isort)
          l=lval_pb(i,isort)
          do ir=0,nrad(isort)
            work(ir)=r(ir,isort)**(l+2)*fun_pb(ir+mtpb)*dr(ir,isort)
          enddo
          mi0(i,isort)=dqdall(h(isort),work,nrad(isort))
        enddo
      enddo			!! over isort
      deallocate(pp,pp1,qq,qq1,work)
      end
