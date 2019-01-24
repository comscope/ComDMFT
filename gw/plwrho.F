      SUBROUTINE PLWRHO(radius,rofi,roatom)
**************************************************************
* Calculates pseudocharge density by fast fourier transform. *
* Written by S.Yu.Savrasov (Lebedev Physical Institute)      *
* Adapted for APW by A.L.Kutepov (VNIITF)                    *
**************************************************************
	use atom_mod
	use manager_mod
      use parallel_mod
      use sitegen_mod
	use solid_mod
	use units_mod
      IMPLICIT none
#ifdef MPI
      include 'mpif.h'
#endif
	real*8, intent(in) :: rofi(0:nrmax,nsort),radius,
     &                      roatom(0:nrmax,nsort)
      INTEGER :: i,isort,id,i0,mfft3,mfft2,mfft1,ifft,iatom,isite,
     &	         iatomr,keyout,nvcmax,nsite
      logical :: jimag
      real*8 :: rhot(5),sqpi4,step,rv,rhos,rhod,rsx,rsy,rsz,rmax,rtx,
     &	        rty,rtz,rt,rtp,rhopnt
	integer, allocatable :: iref(:)
	real*8, allocatable :: ad(:,:),expf(:),pref(:),rhors(:),rhois(:),
     &	                     taut(:,:),site(:,:)
	nvcmax=30000
	allocate(iref(nvcmax),ad(3,3),expf(nsort),pref(nsort),rhors(mmm),
     &	       rhois(mmm),taut(3,natom),site(3,nvcmax))
      allocate(lastg_site(nvcmax/2))
      allocate(iref_vgen(nvcmax))
	allocate(vect_site(3,0:nvcmax))
	allocate(rco_site(nvcmax/2))
	allocate(tauc_site(3,natom))
	allocate(vec_vgen(3,nvcmax))
	allocate(rw_vgen(nvcmax))
      sqpi4=sqrt(4.d0*pi)
      do i=1,3
         ad(i,1)=rbas(i,1)/dfloat(mdiv(1))
         ad(i,2)=rbas(i,2)/dfloat(mdiv(2))
         ad(i,3)=rbas(i,3)/dfloat(mdiv(3))
      enddo
      DO ISORT=1,NSORT
	  if(z(isort)>0.1d0) then
C  DEFINE SMOOTH MT-PART BY A*EXP(-B*R^2)
          STEP=SMT(ISORT)/100.D0
          do id=1,5
            rv=smt(isort)+step*(id-3)
            i0=0
            do i=nr(isort),0,-1
              if(rofi(i,isort).lt.rv) then
                i0=i-1
                exit
              endif
            enddo
            call fint(roatom(i0,isort),rofi(i0,isort),4,rv,rhot(id))
          enddo
          RHOS=RHOT(3)
          RHOD=(RHOT(1)-RHOT(5)+8.D0*(RHOT(4)-RHOT(2)))/12.D0/STEP
          EXPF(ISORT)=-RHOD/RHOS/2.D0/SMT(ISORT)
          PREF(ISORT)=RHOS*EXP(EXPF(ISORT)*SMT(ISORT)**2)
	  else
	    expf(isort)=1.d0
	    pref(isort)=0.d0
	  endif
      enddo
C  DO GLOBALLY OVER  GRID POINTS OF FFT
      rhois=0.d0
      if(me==0) then
        DO MFFT3=0,MDIV(3)-1
          DO MFFT2=0,MDIV(2)-1
            DO MFFT1=0,MDIV(1)-1
              IFFT=1+MFFT1+MFFT2*MDIV(1)+MFFT3*MDIV(1)*MDIV(2)
              RHORS(IFFT)=0.D0
C  RADIUS-VECTOR IS GIVEN BY
              RSX=MFFT1*AD(1,1)+MFFT2*AD(1,2)+MFFT3*AD(1,3)
              RSY=MFFT1*AD(2,1)+MFFT2*AD(2,2)+MFFT3*AD(2,3)
              RSZ=MFFT1*AD(3,1)+MFFT2*AD(3,2)+MFFT3*AD(3,3)
C  CALCULATE ENVIRONMENT SITES FOR THE POINT "RS"
              DO IATOM=1,NATOM
                TAUT(1,IATOM)=TAU(1,IATOM)-RSX
                TAUT(2,IATOM)=TAU(2,IATOM)-RSY
                TAUT(3,IATOM)=TAU(3,IATOM)-RSZ
              enddo
              RMAX=RADIUS/PAR+1.D-05
              call sitegen(rmax,rbas(1,1),rbas(1,2),rbas(1,3),taut,
     &                     natom,nsite,site,iref,nvcmax)
              DO isite=1,nsite
                IATOMR=iref(isite)
                ISORT=IS(IATOMR)
                RTX=-site(1,isite)
                RTy=-site(2,isite)
                RTz=-site(3,isite)
                RT=SQRT(RTX*RTX+RTY*RTY+RTZ*RTZ)
                KEYOUT=1
                rtp=rt*par
                IF(rtp.LE.SMT(ISORT)) KEYOUT=0
C  CALCULATE PSEUDOCHARGE DENSITY AT THE POINT "RT"
                IF(KEYOUT.EQ.1)THEN
                  i0=0
                  do i=nr(isort),0,-1
                    if(rofi(i,isort).lt.rtp) then
                      i0=i-1
                      goto 6
                    endif
                  enddo
6                 call fint(roatom(i0,isort),rofi(i0,isort),4,rtp,
     &                      rhopnt)
                ELSEIF(KEYOUT.EQ.0)THEN
                  RHOPNT=PREF(ISORT)*EXP(-EXPF(ISORT)*rtp**2)
                ENDIF
                RHORS(IFFT)=RHORS(IFFT)+RHOPNT/sqpi4
              enddo
            enddo
          enddo
        enddo
      endif
#ifdef MPI
	if(goparr) call brdcst(msgdbl,rhors,8*mmm,master,MPI_COMM_WORLD)
#endif	
      call timel('******* fourier started ************')
      if(complex_ro) then
        jimag=.true.
        call fourier(rhors,rhois,indstar,rointr,0,igbsf,gbs,0,0,
     &               .true.,jimag)
        call symscal(rointr)
      else
        jimag=.false.
        call fourier(rhors,rhois,indstar,rointr,0,igbsf,gbs,0,0,
     &               .true.,jimag)
        call symscal(rointr)
      endif
	deallocate(iref,ad,expf,pref,rhors,rhois,taut,site)
      deallocate(lastg_site,vect_site,rco_site,tauc_site,iref_vgen,
     &           vec_vgen,rw_vgen)
      end
