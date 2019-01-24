      subroutine setfgtl
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: i0,istar,i,i1,i2,isort,iatom,l,m,lm,ind,irad,
     &             ind_istar
      real*8 :: pi2,pi4,sqpi4,g,fas,const,ss,amegamt,rv,con,cg0,cs,ddot
      complex*16 :: c1,c3 
	real*8, allocatable :: ylm(:,:),bess(:,:),bessd(:),work(:),
     &	                     work1(:),work2(:)
	allocate(ylm(limlpb,maxstar),bess(0:maxpb+1,0:maxnrad),
     &	       bessd(0:maxpb+1),work(0:maxnrad),work1(0:maxnrad),
     &         work2(0:maxnrad))
      pi2=pi+pi
      pi4=pi2+pi2
      sqpi4=sqrt(pi4)
      i0=indstar(n_mpi_istar(me+1))
      do ind_istar=1,ndim_istar(me+1)
        istar=n_mpi_istar(me+1)+ind_istar
        g=rcostar(istar)
        if(istar.gt.1) then
          do i=indstar(istar-1)+1,indstar(istar)
            i1=i-indstar(istar-1)
            call sphharm(gbs(1,i)/g,gbs(2,i)/g,gbs(3,i)/g,maxpb,
     &                   ylm(1,i1))
          enddo
        endif
        do isort=1,nsort
          do iatom=1,natom
            if(is(iatom).eq.isort) goto 11
	    enddo
11	    continue
          do l=0,lmpb(isort)
            c3=pi4*dcmplx(0.d0,1.d0)**l
            do m=-l,l
              lm=l*(l+1)+m+1
              ind=indfgtl(lm,isort)
              if(ind.eq.0) goto 1
              if(istar.eq.1) then
                if(l.eq.0) then
                  agtlr(1,isort,1)=sqpi4
                else
                  agtlr(ind,isort,1)=0.d0
                endif
                agtli(ind,isort,1)=0.d0
              else
                do i=indstar(istar-1)+1,indstar(istar)
                  i1=i-indstar(istar-1)
                  i2=i-i0
                  fas=pi2*ddot(3,gbs(1,i),1,tau(1,iatom),1)
                  c1=dcmplx(cos(fas),sin(fas))*c3
                  agtlr(ind,isort,i2)=dble(c1)*ylm(lm,i1)
                  agtli(ind,isort,i2)=dimag(c1)*ylm(lm,i1)
                enddo
              endif
1             continue
            enddo   !!! over m
          enddo   !!! over l
        enddo   !!! over isort
      enddo   !!! over istar
	deallocate(ylm)
      do ind_istar=1,ndim_istar(me+1)
        istar=n_mpi_istar(me+1)+ind_istar
        g=rcostar(istar)*pi2/par
        do isort=1,nsort
          do iatom=1,natom
            if(is(iatom).eq.isort) goto 12
	    enddo
12	    continue
          const=1.d0/pi4/amega
          ss=smt(isort)**2
          amegamt=pi4*smt(isort)**3/3.d0
          do irad=0,nrad(isort)
            rv=r(irad,isort)
            call BESSR(g,rv,bess(0,irad),bessd,lmpb(isort)+1)
            rv=rv*rv
            work(irad)=(1.d0-rv/ss)**nfgtl
          enddo
          do l=0,lmpb(isort)
            con=const/smt(isort)**(l+3)
            if(istar.ne.1) cg0=sqpi4*smt(isort)**(l+2)*
     &                           bess(l+1,nrad(isort))/g
            do irad=0,nrad(isort)
              work1(irad)=weight(irad,isort)*work(irad)*bess(l,irad)
              work(irad)=work(irad)*r(irad,isort)/smt(isort)
            enddo
            call dqd(h(isort),work1,work2,nrad(isort)+1)
            do m=-l,l
              lm=l*(l+1)+m+1
              ind=indfgtl(lm,isort)
              if(ind.eq.0) goto 2
              cs=con*work2(nrad(isort))
              do i=indstar(istar-1)+1,indstar(istar)
                i2=i-i0
                fgtlr(ind,isort,i2)=cs*agtlr(ind,isort,i2)
                fgtli(ind,isort,i2)=-cs*agtli(ind,isort,i2)
                if(istar.eq.1) then
                  if(l.eq.0) then
                    bgtlr(1,1,isort)=amegamt
                    cgtlr(1,1,isort)=sqpi4
                  else
                    bgtlr(1,ind,isort)=0.d0
                    cgtlr(1,ind,isort)=0.d0
                  endif
                  bgtli(1,ind,isort)=0.d0
                  cgtli(1,ind,isort)=0.d0
                else
                  bgtlr(i2,ind,isort)=cg0*agtlr(ind,isort,i2)
                  cgtlr(i2,ind,isort)=bess(l,nrad(isort))*
     &                                      agtlr(ind,isort,i2)
                  bgtli(i2,ind,isort)=cg0*agtli(ind,isort,i2)
                  cgtli(i2,ind,isort)=bess(l,nrad(isort))*
     &                                      agtli(ind,isort,i2)
                endif
              enddo   !!! over i
2             continue
            enddo   !!! over m
          enddo   !!! over l
        enddo   !!! over isort
      enddo   !!! over istar
	deallocate(bess,bessd,work,work1,work2)
      end
