      subroutine smultro
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
	integer :: i,isort,l,m,lm,isym,ind,mt,mt1,irad,istar,i2,ndimv,
     &	         iatom,iat,ig,igb1,igb2,igb3,ii2,nn,i0,iplw0,j,
     &           j2,ipr,ind_istar
	real*8 :: sqpi4,pi2,c1,csq0,const,amro,amropw,det,det3,fas,cc,cs,
     &	        rv,csq,fac2,qq(3)
      real*8, allocatable :: work2(:),work1(:),qpsro(:,:),psror(:),
     &                       psroi(:)
	allocate(work2(0:maxnrad),work1(0:maxnrad),qpsro(maxnsym,nsort),
     &	       psror(nplwrop))
      if(complex_ro) allocate(psroi(nplwrop))
c   calculate moments of the full density within sphere
      sqpi4=sqrt(4.d0*pi)
	pi2=pi+pi
      c1=2.d0
      do i=2,nfgtl
         c1=c1*2.d0*i
      enddo
      csq0=dfloat(nspin)/2.d0
      csq=csq0/sqpi4
      iplw0=indstar(n_mpi_istar(me+1))
      do isort=1,nsort
        do l=0,lmpb(isort)
          const=fac2(l+l+nfgtl+nfgtl+3)/c1/fac2(l+l+1)/sqpi4
          do m=-l,l
            lm=l*(l+1)+m+1
            if(sym(lm,isort))then
              isym=lmsym(lm,isort)
              ind=indfgtl(lm,isort)
              mt=indmt(isym,isort,1)
              mt1=indmt(isym,isort,nspin)
              do irad=0,nrad(isort)
                work1(irad)=(ro(mt+irad)+ro(mt1+irad))*csq*
     &                            rl(irad,l,isort)*weight(irad,isort)
              enddo
              call dqd(h(isort),work1,work2,nrad(isort)+1)
              amro=work2(nrad(isort))
              if(lm.eq.1) amro=amro-z(isort)
              amropw=0.d0
              do ind_istar=1,ndim_istar(me+1)
                istar=n_mpi_istar(me+1)+ind_istar
                do i=indstar(istar-1)+1,indstar(istar)
                  i2=i-iplw0
                  amropw=amropw+real(rointr(i,1)+rointr(i,nspin))
     &                         *bgtlr(i2,ind,isort)*csq0
                  if(complex_ro) amropw=amropw
     &                                 -imag(rointr(i,1)+
     &                                       rointr(i,nspin))*csq0*
     &                                  bgtli(i2,ind,isort)
                enddo
              enddo
#ifdef MPI
              call DGOP(amropw,1,'  +',MPI_COMM_WORLD)
#endif              
              qpsro(isym,isort)=(amro-amropw)*const
            endif
          enddo  !!! over m
        enddo  !!! over l
      enddo  !!! over isort
      det=det3(gbas(1,1),gbas(1,2),gbas(1,3))
      do ind_istar=1,ndim_istar(me+1)
        istar=n_mpi_istar(me+1)+ind_istar
        do i=indstar(istar-1)+1,indstar(istar)
          i2=i-iplw0
          psror(i2)=0.d0
          if(complex_ro) psroi(i2)=0.d0
          do isort=1,nsort
            ndimv=(lmpb(isort)+1)**2
c --------- First atoms of the classes ----------------------------------
            iatom=iat_1(isort)
            do iat=1,natom
	        do ig=1,ngroup
	          if(ip(iat,ig).eq.iatom) then
	            fas=pi2*dot_product(gbs(:,i),shift(:,ig))
	            cc=cos(fas)
	            cs=sin(fas)
	            call rotate(gbs(1,i),gbs(2,i),gbs(3,i),qq(1),qq(2),
     &				          qq(3),u(2,ig),1)    !!  qq=R^(-1) G
                  igb1=idnint(det3(qq,gbas(1,2),gbas(1,3))/det)
                  igb2=idnint(det3(gbas(1,1),qq,gbas(1,3))/det)
                  igb3=idnint(det3(gbas(1,1),gbas(1,2),qq)/det)
	            ii2=indplw(igb1,igb2,igb3)-iplw0
                  do lm=1,ndimv
                    if(sym(lm,isort))then
                      isym=lmsym(lm,isort)
                      ind=indfgtl(lm,isort)
                      psror(i2)=psror(i2)+qpsro(isym,isort)
     &			             *(cc*fgtlr(ind,isort,ii2)
     &                          +cs*fgtli(ind,isort,ii2))
                      if(complex_ro) psroi(i2)=psroi(i2)
     &                                        +qpsro(isym,isort)
     &			                              *(cc*fgtli(ind,isort,ii2)
     &                                         -cs*fgtlr(ind,isort,ii2))
                    endif
                  enddo
	            goto 2
	          endif
	        enddo   !! over ig
2             continue
	      enddo    !! over iat of the given sort
          enddo      !! over isort
        enddo   !! over G of the given star
      enddo     !! over stars
      const=2.d0*par*par/pi
      v_intr_h(1)=(0.d0,0.d0)
      do ind_istar=1,ndim_istar(me+1)
        istar=n_mpi_istar(me+1)+ind_istar
        if(istar.ne.1) then
          do j=indstar(istar-1)+1,indstar(istar)
            j2=j-iplw0
            rv=psror(j2)
            rv=rv+real(rointr(j,1)+rointr(j,nspin))*csq0
            v_intr_h(j)=const*rv/rcostar(istar)**2
            if(complex_ro) then
              rv=psroi(j2)
              rv=rv+imag(rointr(j,1)+rointr(j,nspin))*csq0
              v_intr_h(j)=v_intr_h(j)+(0.d0,1.d0)*const*rv
     &                                /rcostar(istar)**2
            endif
          enddo
        endif
      enddo
      if(goparr) then
        do ipr=0,nproc-1
          if(ipr.eq.me) then
            nn=nplwrop
            i0=iplw0
          endif
#ifdef MPI
          call brdcst(msgint,nn,4,ipr,MPI_COMM_WORLD)
#endif          
          if(nn==0) exit
#ifdef MPI
          call brdcst(msgint,i0,4,ipr,MPI_COMM_WORLD)
          call brdcst(msgdbl,v_intr_h(i0+1),nn*16,ipr,MPI_COMM_WORLD)
#endif          
        enddo
      endif
      vatmt = 0.0d0
      do isort=1,nsort
        ndimv=(lmpb(isort)+1)**2
        do lm=1,ndimv
          if(sym(lm,isort))then
            isym=lmsym(lm,isort)
            ind=indfgtl(lm,isort)
            vatmt(isym,isort)=0.d0
            do ind_istar=1,ndim_istar(me+1)
              istar=n_mpi_istar(me+1)+ind_istar
              if(istar==1) cycle
              do i=indstar(istar-1)+1,indstar(istar)
                i2=i-iplw0
                vatmt(isym,isort)=vatmt(isym,isort)
     &                           +real(v_intr_h(i))
     &                           *cgtlr(i2,ind,isort)
              enddo
            enddo
            if(complex_ro) then
              do ind_istar=1,ndim_istar(me+1)
                istar=n_mpi_istar(me+1)+ind_istar
                if(istar==1) cycle
                do i=indstar(istar-1)+1,indstar(istar)
                  i2=i-iplw0
                  vatmt(isym,isort)=vatmt(isym,isort)
     &                             -imag(v_intr_h(i))
     &                             *cgtli(i2,ind,isort)
                enddo
              enddo
            endif
          endif
        enddo   !!! over lm
      enddo      !!! over isort
#ifdef MPI
      call DGOP(vatmt,maxnsym*nsort,'  +',MPI_COMM_WORLD)
#endif      
	deallocate(work2,work1,qpsro,psror)
      if(complex_ro) deallocate(psroi)
      call timel('**** smultro finished **************')
1000  format(1x,'Vmt0 = ',e15.8)
1010  format(1x,'After ',i6,' stars Vmt0 = ',e15.8)
1020  format(1x,'Isort Isym V_at_MT Psatmt = ',2i4,2e15.8)
1030  format(1x,'Isort Isym AMRO AMROPW = ',2i4,2e15.8)
1040  format(1x,'Rointr Psror Vintr = ',i4,3e15.8)
1050  format(1x,'Istar Qpsro Fgtlr Psror = ',i4,3e13.6)
      end
