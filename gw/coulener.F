      subroutine coulener(rovhartr)
**************************************************************
*  Calculate electrostatic contribution to the total energy  *
**************************************************************
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
	real*8, intent(out) :: rovhartr
      integer :: isort,isym,mt,mt1,irad,ispin,i,j,ia,ib,ic,ind,l,m,lm,
     &           ishplw
      real*8 :: ecouint,sqpi4,const,ecoumt,tr1,tr2,sum,b,dqdall
	real*8, allocatable :: work1(:)
	allocate(work1(0:maxnrad))
      sqpi4=2.d0*sqrt(pi)
      ishplw=nplwro/nproc
      if(ishplw*nproc.lt.nplwro) ishplw=ishplw+1
#ifdef MPI
	if(goparr) then
	  call brdcst(msgdbl,ro_atom,8*(nrmax+1)*nsort,master,
     &	            MPI_COMM_WORLD)
	  call brdcst(msgdbl,eny,8*maxntle*(nrel*maxb+1)*nsort*nspin,
     &	            master,MPI_COMM_WORLD)
	  call brdcst(msgdbl,e_core,8*ncormax*nsort*nspin,
     &	            master,MPI_COMM_WORLD)
	  call brdcst(msgdbl,ro_core,8*(nrmax+1)*nsort*nspin,
     &	            master,MPI_COMM_WORLD)
      endif
#endif      
c   calculate mt-contribution to the coulomb energy
      ecoumt=0.d0
      do isort=1,nsort
        do l=0,lmpb(isort)
          do m=-l,l
            lm=l*(l+1)+m+1
            if(sym(lm,isort))then
              isym=lmsym(lm,isort)
              if(l.eq.0)then
                mt=indmt(1,isort,1)
                mt1=indmt(1,isort,nspin)
                const=0.25d0*dfloat(nspin)
                do irad=0,nrad(isort)
                  tr1=v_mt_h(mt+irad)*r(irad,isort)**2*dr(irad,isort)
                  tr2=(ro(mt+irad)+ro(mt1+irad))
                  work1(irad)=tr1*tr2
                enddo   !!! over irad
                sum=const*dqdall(h(isort),work1,nrad(isort))
                const=-sqpi4*z(isort)*dfloat(nspin)/2.d0
                do irad=0,nrad(isort)
                  tr1=r(irad,isort)*dr(irad,isort)
                  tr2=ro(mt+irad)+ro(mt1+irad)
                  work1(irad)=tr1*tr2
                enddo
                sum=sum+const*dqdall(h(isort),work1,nrad(isort))
	          ecoumt=ecoumt+sum*dfloat(nhsort(isort))
              else
                mt=indmt(isym,isort,1)
                mt1=indmt(isym,isort,nspin)
                const=0.25d0*dfloat(nspin)
                do irad=0,nrad(isort)
                  work1(irad)=(ro(mt+irad)+ro(mt1+irad))*v_mt_h(mt+irad)
     &				       *r(irad,isort)**2*dr(irad,isort)
                enddo
                sum=const*dqdall(h(isort),work1,nrad(isort))
                ecoumt=ecoumt+sum*dfloat(nhsort(isort))
              endif
            endif
          enddo   !!! over m
        enddo   !!! over l
      enddo   !!! over isort
      ecouint=0.d0
c  calculate contribution from the interstitial region
        do ispin=1,nspin
          b=0.d0
          do i=1+ishplw*me,min(nplwro,ishplw*(me+1))
            do j=1,nplwro
              ia=igbs(1,i)+igbs(1,j)
              ib=igbs(2,i)+igbs(2,j)
              ic=igbs(3,i)+igbs(3,j)
              ind=indplw(ia,ib,ic)
              b=b+real(rointr(i,ispin))*real(v_intr_h(j))*sovr(ind)
            enddo
          enddo
          if(complex_ro) then
            do i=1+ishplw*me,min(nplwro,ishplw*(me+1))
              do j=1,nplwro
                ia=igbs(1,i)+igbs(1,j)
                ib=igbs(2,i)+igbs(2,j)
                ic=igbs(3,i)+igbs(3,j)
                ind=indplw(ia,ib,ic)
                b=b-imag(rointr(i,ispin))*imag(v_intr_h(j))*sovr(ind)
                b=b-(imag(rointr(i,ispin))*real(v_intr_h(j))+
     &               real(rointr(i,ispin))*imag(v_intr_h(j)))*sovi(ind)
              enddo
            enddo
          endif
#ifdef MPI
          call DGOP(b,1,'  +',MPI_COMM_WORLD)
#endif          
          ecouint=ecouint+0.5d0*b*amega
        enddo
        if(maswrk) write(iun,1010)ecoumt,ecouint
	deallocate(work1)
	rovhartr=2.d0*(ecoumt+ecouint)
1010  format(1x,'Coulomb Energy : mt,int-parts are',2(1x,g17.10))
1011  format(1x,'Coulomb Energy : one,two-site parts are',2(1x,g17.10))
      end