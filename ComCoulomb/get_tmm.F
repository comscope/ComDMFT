      subroutine get_tmm
      

      ! evaluate R=0 component of bare coulomb interaction

      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod

      implicit none
      integer :: isort,l,iej,j,mtpb,iei,i,mtpbi,ir,li,je1,le1,ispin,li1,
     &  l1,je3,le3,l2,le,le2,je,je2
      double precision :: pi8,pi8l,dqdall,deg,fc,c2,ss
      double precision, allocatable ::
     $  work(:),pp(:),tmp(:),f(:),g(:),qq(:),
     &  pp1(:)

      if(obietto=='   solid') then
        allocate(work(0:maxnrad),pp(0:maxnrad),tmp(0:maxnrad),
     &    f(0:maxnrad),g(0:maxnrad))
      else if(obietto=='    atom') then
        allocate(work(0:nrmax),pp(0:nrmax),tmp(0:nrmax),f(0:nrmax),
     &    g(0:nrmax))
      endif

      c2=clight**2
      pi8=8.d0*pi
      tmm_pb=0.d0

      do isort=1,nsort
c     ------ Getting MM-values  ----------------------------------------
        do l=0,lmpb(isort)
          pi8l=pi8/dfloat(l+l+1)
          do iej=1,ntle_pb(l,isort)
            if(obietto=='   solid') then
              j=ind_prod0(iej,l,isort)
              mtpb=ind_pb(j,isort)
              call var_int1(dr(0,isort),fun_pb(mtpb),f,g,tmp(0),
     &          r(0,isort),rl(0,l,isort),work,nrad(isort)+1,
     &          h(isort))
            else if(obietto=='    atom') then
              call var_int1(dr_atom(0),pb_atom_vv(0,iej,l),f,g,tmp(0),
     &          r_atom(0,1),rl_atom(0,l),work,nrmax+1,h(1))
            endif
            do iei=1,ntle_pb(l,isort)
              if(obietto=='   solid') then
                i=ind_prod0(iei,l,isort)
                mtpbi=ind_pb(i,isort)
                do ir=0,nrad(isort)
                  pp(ir)=tmp(ir)*fun_pb(mtpbi+ir)*dr(ir,isort)
     &              *r(ir,isort)**2
                enddo
                tmm_pb(iei,iej,l,isort)=pi8l
     &            *dqdall(h(isort),pp,nrad(isort))
              else if(obietto=='    atom') then
                do ir=0,nrmax
                  pp(ir)=tmp(ir)*pb_atom_vv(ir,iei,l)*dr_atom(ir)
     &              *r_atom(ir,1)**2
                enddo
                tmm_pb(iei,iej,l,1)=pi8l*dqdall(h(1),pp,nrmax)
              endif
            enddo
          enddo
          do iej=1,ntle_pb(l,isort)
            do iei=1,ntle_pb(l,isort)
              tmm_pb(iei,iej,l,isort)=0.5d0*(tmm_pb(iei,iej,l,isort)+
     &          tmm_pb(iej,iei,l,isort))
              tmm_pb(iej,iei,l,isort)=tmm_pb(iei,iej,l,isort)
            enddo
          enddo
        enddo
      enddo                     !!! over isort
      deallocate(work,pp,tmp,f,g)
c     ---------------------------------------------------------
      if(obietto=='    atom') then
        allocate(f(0:nrmax),g(0:nrmax),work(0:nrmax),tmp(0:nrmax),
     &    pp(0:nrmax),qq(0:nrmax))
        allocate(pp1(0:nrmax))
        v4_atom=0.d0
        do ispin=1,nspin
          do li=1,nrel*lmb(1)+1
            if(irel.le.1) then
              l=li-1
              deg=l+l+1.d0
            else
              l=li/2
              deg=dfloat(((li+1)/2)*2)
            endif
            fc=pi8/deg
            do je1=1,ntle(l,1)
              le1=indfun0(1,je1,li,1)
              do ir=0,nrmax
                pp(ir)=pval_atom(ir,le1,ispin)
                if(irel.ge.1) qq(ir)=qval_atom(ir,le1,ispin)
              enddo
              do li1=1,nrel*lmb(1)+1
                if(irel.le.1) then
                  l1=li1-1
                else
                  l1=li1/2
                endif
                do je3=1,ntle(l1,1)
                  le3=indfun0(1,je3,li1,1)
                  do l2=0,maxpb
                    call var_int(dr_atom(0),pval_atom(0,le3,ispin),
     &                qval_atom(0,le3,ispin),pp(0),qq(0),f,g,
     &                tmp(0),r_atom(0,1),rl_atom(0,l2),work,
     &                nrmax+1,h(1))
                    do je=1,ntle(l,1)
                      le=indfun0(1,je,li,1)
                      do je2=1,ntle(l1,1)
                        le2=indfun0(1,je2,li1,1)
                        do ir=0,nrmax 
                          pp1(ir)=pval_atom(ir,le,ispin)
     &                      *pval_atom(ir,le2,ispin)
                          if(irel.ge.1) pp1(ir)=pp1(ir)
     &                      +qval_atom(ir,le,ispin)
     &                      *qval_atom(ir,le2,ispin)/c2
                          pp1(ir)=pp1(ir)*tmp(ir)*dr_atom(ir)
                        enddo
                        ss=dqdall(h(1),pp1,nrmax)*d_lll(l2,li,li1)**2
                        i=ntle(l1,1)*(je2-1)+je
                        j=ntle(l1,1)*(je3-1)+je1
                        v4_atom(i,j,li,li1,ispin)=
     &                    v4_atom(i,j,li,li1,ispin)+fc*ss
                      enddo     !! over je2
                    enddo       !! over je
                  enddo         !! over l2
                enddo           !! je3
              enddo             !! over li1
            enddo               !! over je1
          enddo                 !! over li
        enddo                   !! over ispin
        deallocate(f,g,work,tmp,pp,qq,pp1)
      endif
      end
