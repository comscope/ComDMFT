      subroutine matchav(f,fintr,jreal,jimag)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      logical :: jreal,jimag
      real*8, intent(in) :: f(maxmtb)
      complex*16, intent(in) :: fintr(3,nplwro)
	integer :: iplw0,i,iend,isort,istar,l,lm,m,isym,ind,mt,j,j2,
     &             ind_istar
	real*8 :: const,rco,all,valmt,dval1t,dval1mt,dval2t,dval2mt,
     &	      valint,dval1int,dval2int,ss,e
	real*8, allocatable :: bess(:,:),bessd(:,:),bessd2(:,:)
	allocate(bess(0:maxpb,nstar),bessd(0:maxpb,nstar),
     &	     bessd2(0:maxpb,nstar))
      const=2.d0*pi/par
      iplw0=indstar(n_mpi_istar(me+1))
      iend=0
      do i=1,72
        if(text(i:i).eq.'!') then
          iend=i
          goto 1
        endif
      enddo
1     if(maswrk) write(iun,'(72a1)')(text(i:i),i=1,iend)
      do isort=1,nsort
        if(maswrk) write(iun,'(2a10,i4)')txtel(isort),'  Isort = ',isort
        do istar=1,nstar
          rco=const*rcostar(istar)
          call bessr(rco,smt(isort),bess(0,istar),bessd(0,istar),
     &               lmpb(isort))
          e=rco**2
          do l=0,lmpb(isort)
            all=dfloat(l*(l+1))/smt(isort)**2
            bessd2(l,istar)=(all-e)*bess(l,istar)-2.d0*bessd(l,istar)
     &                                                    /smt(isort)
          enddo
        enddo
        do i=1,3
          do l=0,lmpb(isort)
            do m=-l,l
              lm=l*(l+1)+m+1
              if(symb(lm,i,isort)) then
                isym=lmsymb(lm,i,isort)
                ind=indfgtl(lm,isort)
                mt=indmtb(isym,isort)+nrad(isort)
                valmt=f(mt)
                dval1t=(11.d0*f(mt)-18.d0*f(mt-1)+9.d0*f(mt-2)-
     &                                  2.d0*f(mt-3))/6.d0/h(isort)
                dval1mt=dval1t/dr(nrad(isort),isort)
                dval2t=(2.d0*f(mt)-5.d0*f(mt-1)+4.d0*f(mt-2)-
     &                                  f(mt-3))/h(isort)**2
                dval2mt=(dval2t-dval1t)/dr(nrad(isort),isort)**2
                valint=0.d0
                dval1int=0.d0
                dval2int=0.d0
                do ind_istar=1,ndim_istar(me+1)
                  istar=n_mpi_istar(me+1)+ind_istar
                  ss=0.d0
                  do j=indstar(istar-1)+1,indstar(istar)
                    j2=j-iplw0
                    if(jreal) ss=ss+real(fintr(i,j))*agtlr(ind,isort,j2)
                    if(jimag) ss=ss-imag(fintr(i,j))*agtli(ind,isort,j2)
                  enddo
                  valint=valint+ss*bess(l,istar)
                  dval1int=dval1int+ss*bessd(l,istar)
                  dval2int=dval2int+ss*bessd2(l,istar)
                enddo
#ifdef MPI
                call DGOP(valint,1,'  +',MPI_COMM_WORLD)
                call DGOP(dval1int,1,'  +',MPI_COMM_WORLD)
                call DGOP(dval2int,1,'  +',MPI_COMM_WORLD)
#endif                
                if(maswrk) then
                  if(i.eq.1) write(iun,*)' For X : L = ',l,' M = ',m
                  if(i.eq.2) write(iun,*)' For Y : L = ',l,' M = ',m
                  if(i.eq.3) write(iun,*)' For Z : L = ',l,' M = ',m
                  write(iun,1000)valmt,valint
                  write(iun,1001)dval1mt,dval1int
                  write(iun,1002)dval2mt,dval2int
                endif
              endif
            enddo   !!! over m
          enddo   !!! over l
        enddo   !!! over i
      enddo   !!! over isort
1000  format(' Value    Mt : ',f15.10,'  Interst. : ',f15.10)
1001  format(' d / dr   Mt : ',f15.10,'  Interst. : ',f15.10)
1002  format(' d2/ dr2  Mt : ',f15.10,'  Interst. : ',f15.10)
	deallocate(bess,bessd,bessd2)
      end
