      subroutine get_cfgt_old(qg2,qg,ig,vmi)
	use atom_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ig
      real*8, intent(in) :: qg2,qg(3)
      complex*16, intent(out) :: vmi(n_pbmt,*)
	integer :: iatom,isort,irad,lm,l,m,i0,ie,j,mtpb,km
      real*8 :: pi4,pi2,sqpi4,qg1,fas,rv,cs,ddot,dqdall
	complex*16 c1,c3,agt
	real*8, allocatable :: ylm(:),bess(:,:),bessd(:),work(:),work1(:)
	allocate (ylm(limlpb),bess(0:maxpb,0:maxnrad),
     &	      bessd(0:maxpb),work(0:maxnrad),work1(0:maxnrad))
	pi4=4.d0*pi
	pi2=pi+pi
	sqpi4=sqrt(pi4)
	qg1=sqrt(qg2)
	if(qg1.gt.1.d-6) then
	  call sphharm(qg(1)/qg1,qg(2)/qg1,qg(3)/qg1,maxpb,ylm)
	else
	  ylm=0.d0
	  ylm(1)=1.d0/sqpi4
	endif
	do iatom=1,natom
	  i0=iopb(iatom)-1
	  isort=is(iatom)
        fas=par*ddot(3,qg,1,tau(1,iatom),1)
        c1=exp(dcmplx(0.d0,fas))
        do irad=0,nrad(isort)
          rv=r(irad,isort)
          call BESSR(qg1,rv,bess(0,irad),bessd,lmpb(isort))
          rv=rv*rv
          work(irad)=rv*dr(irad,isort)
        enddo
        do l=0,lmpb(isort)
          if(qg1.gt.1.d-6) then
            c3=pi4*dcmplx(0.d0,1.d0)**l*c1
	    else
	      c3=dcmplx(pi4,0.d0)
	    endif
	    do ie=1,ntle_pb(l,isort)
	      j=ind_prod0(ie,l,isort)
	      mtpb=ind_pb(j,isort)
            do irad=0,nrad(isort)
              work1(irad)=work(irad)*bess(l,irad)*fun_pb(mtpb+irad)
            enddo
            cs=dqdall(h(isort),work1,nrad(isort))
            agt=c3*cs
            do m=-l,l
	        lm=l*(l+1)+m+1
	        km=ind_prod(ie,lm,isort)
              vmi(i0+km,ig)=ylm(lm)*agt
	      enddo
	    enddo
	  enddo  !! over l
	enddo !! over iatom
	deallocate (ylm,bess,bessd,work,work1)
      end