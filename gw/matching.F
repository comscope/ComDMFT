      subroutine matching(key,f,fintr)
c     key=0   f1
c     key=1   0.5*(f1-f2)
c     key=2   f1+f2
c     key=3   f1-f2
c     key=4   0.5*(f1+f2)
c     key=5   f2
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	implicit none
#ifdef MPI
      include 'mpif.h'
#endif
	integer, intent(in) :: key
	real*8, intent(in) :: f(*)
	complex*16, intent(in) :: fintr(nplwro,*)
      integer :: i,iplw0,ima,icorr,iend,isort,istar,l,m,isym,ind,mt,j,
     &	       j2,mt1,lm,ind_istar
      real*8 :: const,sqpi4,con,rco,e,all,valmt,dval1t,dval1mt,dval2t,
     &	      dval2mt,valint,dval1int,dval2int,ss,vlmt,dvl1mt,dvl2mt,
     &          vlint,dvl1int,dvl2int
	real*8, allocatable :: bess(:,:),bessd(:,:),bessd2(:,:)
	allocate(bess(0:maxpb,nstar),bessd(0:maxpb,nstar),
     &	     bessd2(0:maxpb,nstar))
      const=2.d0*pi/par
      sqpi4=sqrt(4.d0*pi)
      con=2.d0*sqpi4
      iplw0=indstar(n_mpi_istar(me+1))
      ima=1
      if(key.eq.5) ima=2
      icorr=0
      do i=1,60
         if(text(i:i+5).eq.'Coulon') icorr=1
         if(text(i:i+5).eq.'Full P') icorr=1
      enddo
      iend=0
      do i=1,72
         if(text(i:i).eq.'!') then
            iend=i
            goto 1
         endif
      enddo
1     if(maswrk) write(iun,'(72a1)')(text(i:i),i=1,iend)
      do isort=1,nsort
        if(maswrk)write(iun,'(a4,a9,i4)')txtel(isort),' Isort = ',isort
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
        do l=0,lmpb(isort)
          do m=-l,l
            lm=l*(l+1)+m+1
            if(sym(lm,isort)) then
              isym=lmsym(lm,isort)
              ind=indfgtl(lm,isort)
              mt=indmt(isym,isort,ima)+nrad(isort)
              valmt=f(mt)
              if(l.eq.0.and.icorr.eq.1) valmt=valmt-con*z(isort)
     &			                                 /smt(isort)
              dval1t=(11.d0*f(mt)-18.d0*f(mt-1)+9.d0*f(mt-2)-
     &                               2.d0*f(mt-3))/6.d0/h(isort)
              dval1mt=dval1t/dr(nrad(isort),isort)
              if(l.eq.0.and.icorr.eq.1) dval1mt=dval1mt+
     &                                    con*z(isort)/smt(isort)**2
              dval2t=(2.d0*f(mt)-5.d0*f(mt-1)+4.d0*f(mt-2)-
     &                               f(mt-3))/h(isort)**2
              dval2mt=(dval2t-dval1t)/dr(nrad(isort),isort)**2
              if(l.eq.0.and.icorr.eq.1) dval2mt=dval2mt-2.d0*
     &                                    con*z(isort)/smt(isort)**3
              valint=0.d0
              dval1int=0.d0
              dval2int=0.d0
              do ind_istar=1,ndim_istar(me+1)
                istar=n_mpi_istar(me+1)+ind_istar
                ss=0.d0
                do j=indstar(istar-1)+1,indstar(istar)
                  j2=j-iplw0
                  ss=ss+real(fintr(j,ima))*agtlr(ind,isort,j2)
                  if(complex_ro) ss=ss-imag(fintr(j,ima))
     &                                *agtli(ind,isort,j2)
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
              if(key.gt.0.and.key.lt.5) then
                mt1=indmt(isym,isort,2)+nrad(isort)
                vlmt=f(mt1)
                if(l.eq.0.and.icorr.eq.1) vlmt=vlmt-
     &                                    con*z(isort)/smt(isort)
                dval1t=(11.d0*f(mt1)-18.d0*f(mt1-1)+9.d0*f(mt1-2)-
     &                               2.d0*f(mt1-3))/6.d0/h(isort)
                dvl1mt=dval1t/dr(nrad(isort),isort)
                if(l.eq.0.and.icorr.eq.1) dvl1mt=dvl1mt+
     &                                    con*z(isort)/smt(isort)**2
                dval2t=(2.d0*f(mt1)-5.d0*f(mt1-1)+4.d0*f(mt1-2)-
     &                               f(mt1-3))/h(isort)**2
                dvl2mt=(dval2t-dval1t)/dr(nrad(isort),isort)**2
                if(l.eq.0.and.icorr.eq.1) dvl2mt=dvl2mt-2.d0*
     &                                    con*z(isort)/smt(isort)**3
                vlint=0.d0
                dvl1int=0.d0
                dvl2int=0.d0
                do ind_istar=1,ndim_istar(me+1)
                  istar=n_mpi_istar(me+1)+ind_istar
                  ss=0.d0
                  do j=indstar(istar-1)+1,indstar(istar)
                    j2=j-iplw0
                    ss=ss+real(fintr(j,2))*agtlr(ind,isort,j2)
                    if(complex_ro) ss=ss-imag(fintr(j,2))
     &                                  *agtli(ind,isort,j2)
                  enddo
                  vlint=vlint+ss*bess(l,istar)
                  dvl1int=dvl1int+ss*bessd(l,istar)
                  dvl2int=dvl2int+ss*bessd2(l,istar)
	          enddo
#ifdef MPI
                call DGOP(vlint,1,'  +',MPI_COMM_WORLD)
                call DGOP(dvl1int,1,'  +',MPI_COMM_WORLD)
                call DGOP(dvl2int,1,'  +',MPI_COMM_WORLD)
#endif                
                if(key.eq.1) then
                  valmt=0.5d0*(valmt-vlmt)
                  dval1mt=0.5d0*(dval1mt-dvl1mt)
                  dval2mt=0.5d0*(dval2mt-dvl2mt)
                  valint=0.5d0*(valint-vlint)
                  dval1int=0.5d0*(dval1int-dvl1int)
                  dval2int=0.5d0*(dval2int-dvl2int)
                else if(key.eq.2) then
                  valmt=valmt+vlmt
                  dval1mt=dval1mt+dvl1mt
                  dval2mt=dval2mt+dvl2mt
                  valint=valint+vlint
                  dval1int=dval1int+dvl1int
                  dval2int=dval2int+dvl2int
                else if(key.eq.3) then
                  valmt=valmt-vlmt
                  dval1mt=dval1mt-dvl1mt
                  dval2mt=dval2mt-dvl2mt
                  valint=valint-vlint
                  dval1int=dval1int-dvl1int
                  dval2int=dval2int-dvl2int
                else if(key.eq.4) then
                  valmt=0.5d0*(valmt+vlmt)
                  dval1mt=0.5d0*(dval1mt+dvl1mt)
                  dval2mt=0.5d0*(dval2mt+dvl2mt)
                  valint=0.5d0*(valint+vlint)
                  dval1int=0.5d0*(dval1int+dvl1int)
                  dval2int=0.5d0*(dval2int+dvl2int)
                endif
              endif
              if(maswrk) then
                write(iun,*)' For L = ',l,' M = ',m
                write(iun,1000)valmt,valint
                write(iun,1001)dval1mt,dval1int
                write(iun,1002)dval2mt,dval2int
              endif
	      endif
          enddo   !!! over m
        enddo   !!! over l
      enddo   !!! over isort
1000  format(' Value    Mt : ',f15.10,'  Interst. : ',f15.10)
1001  format(' d / dr   Mt : ',f15.10,'  Interst. : ',f15.10)
1002  format(' d2/ dr2  Mt : ',f15.10,'  Interst. : ',f15.10)
	deallocate(bess,bessd,bessd2)
      end
