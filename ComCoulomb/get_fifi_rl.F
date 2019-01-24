      subroutine get_fifi_rl
      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer :: ispin,isort,le1,mt1,mt,le,ir,l
      double precision :: c2,dqdall
      double precision, allocatable :: pp(:),pp1(:),qq(:),qq1(:),work(:)
      allocate(pp(0:maxnrad),pp1(0:maxnrad),qq(0:maxnrad),
     $  qq1(0:maxnrad),work(0:maxnrad))
      c2=clight**2

      do ispin=1,nspin
        do isort=1,nsort
          do le1=1,lfun(isort)
            mt1=ind_wf(le1,isort)
            do ir=0,nrad(isort)
              pp1(ir)=gfun(mt1+ir,ispin)*r(ir,isort)
              if(irel.ge.1) qq1(ir)=gfund(mt1+ir,ispin)*r(ir,isort)
            enddo
            do le=1,lfun(isort)
              mt=ind_wf(le,isort)
              do ir=0,nrad(isort)
                pp(ir)=gfun(mt+ir,ispin)*r(ir,isort)
                if(irel.ge.1) qq(ir)=gfund(mt+ir,ispin)*r(ir,isort)
              enddo
              do l=0,lmpb(isort)
                do ir=0,nrad(isort)
                  work(ir)=pp(ir)*pp1(ir)
                  if(irel>=1) work(ir)=work(ir)+qq(ir)*qq1(ir)/c2
                  work(ir)=work(ir)*rl(ir,l,isort)*dr(ir,isort)
                enddo
                fifi_rl(le,le1,l,isort,ispin)=dqdall(h(isort),work,
     $            nrad(isort))
              enddo		!! over l
            enddo		!! over le
          enddo                 !! over le1
        enddo			!! over isort
      enddo			!! over ispin
      deallocate(pp,pp1,qq,qq1,work)
      end
