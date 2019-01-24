      subroutine mq_prepare(mq,q)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      real*8, intent(in) :: q(3)
      complex*16, intent(out) :: mq(n_pbmtm,n_pbmtm,natom)
	integer :: isort,ir,l,mt,mt1,lf,lf1,ndimb,lm1,lm,iatom,
     &	         ii,iii,l1,lget,nm1,nm,indx,icg1,icg2,icg,lmi,l_i,ie1,
     &           km1,ie,km
      real*8 :: pi4,qq,rv,dqdall,c2,sqpi4,pp,phase,pi2,pi2a
      complex*16 :: cc
	real*8, allocatable :: ylm(:),bess(:,:),bessd(:),work(:),work1(:),
     &                       fi1(:,:,:)	
	allocate (ylm((2*maxpb+1)**2),bess(0:2*maxpb,0:maxnrad),
     &	        bessd(0:2*maxpb),work(0:maxnrad),work1(0:maxnrad))
      allocate(fi1(maxlfpb,maxlfpb,0:2*maxpb))
      pi2=pi+pi
	pi2a=pi2/par
	pi4=4.d0*pi
	sqpi4=sqrt(pi4)
	c2=clight**2
	mq=(0.d0,0.d0)
	qq=sqrt(dot_product(q,q))
	if(qq.gt.1.d-9) then
	  call sphharm(q(1)/qq,q(2)/qq,q(3)/qq,2*maxpb,ylm)
	else
	  ylm=0.d0
	  ylm(1)=1.d0/sqpi4
	endif
	do iatom=1,natom
	  isort=is(iatom)
        do ir=0,nrad(isort)
          rv=r(ir,isort)
          call BESSR(qq*pi2a,rv,bess(0,ir),bessd,2*lmpb(isort))
          rv=rv*rv
          work(ir)=rv*dr(ir,isort)
        enddo
        do l=0,2*lmpb(isort)
	    do lf1=1,lfun_pb(isort)
	      mt1=ind_pb(lf1,isort)
	      do lf=1,lfun_pb(isort)
	        mt=ind_pb(lf,isort)
	        do ir=0,nrad(isort)
	          pp=fun_pb(mt+ir)*fun_pb(mt1+ir)
                work1(ir)=work(ir)*bess(l,ir)*pp
	        enddo
              fi1(lf,lf1,l)=pi4*dqdall(h(isort),work1,nrad(isort))
            enddo   !! over lf
          enddo   !! over lf1
	  enddo  !! over l
	  ndimb=(lmpb(isort)+1)**2
	  do lm1=1,ndimb
	    l1=lget(lm1)
	    nm1=ntle_pb(l1,isort)
	    do lm=1,ndimb
	      l=lget(lm)
	      nm=ntle_pb(l,isort)
	      ii = max0(lm1,lm)
	      iii = min0(lm1,lm)
            indx = (ii*(ii-1))/2 + iii
            icg1 = indxcg(indx)
            icg2 = indxcg(indx+1) - 1 
            do icg = icg1, icg2
	        lmi = jcg(icg)
	        l_i=lget(lmi)
	        cc=(0.d0,1.d0)**l_i*ylm(lmi)
	        do ie1=1,nm1
	          km1=ind_prod(ie1,lm1,isort)
	          lf1=ind_prod0(ie1,l1,isort)
		        do ie=1,nm
	            km=ind_prod(ie,lm,isort)
	            lf=ind_prod0(ie,l,isort)
	            mq(km,km1,iatom)=mq(km,km1,iatom)
     &	                          +cc*cg(icg)*fi1(lf,lf1,l_i)
	          enddo   !! over ie
	        enddo   !! over ie1
	      enddo   !! over icg
	    enddo   !! over lm
	  enddo  !! over lm1
	  phase=pi2*dot_product(q,tau(:,iatom))
	  mq(:,:,iatom)=dcmplx(cos(phase),sin(phase))*mq(:,:,iatom)
	enddo !! over iatom
	deallocate (ylm,bess,bessd,work,work1,fi1)
      end