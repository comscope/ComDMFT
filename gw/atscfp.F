      subroutine atscfp(isort,nr0,rofi,roatom,rmax)
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
	integer, intent(in) :: isort,nr0
      real*8,intent(in) :: rofi(0:nr0)
      real*8,intent(out) :: roatom(0:nr0)
	integer :: niter,ir,ic,li,ie,nval,l,ndim,nds
	real*8 :: vhrmax,sqpi4,tol,tolrsq,beta,b,rmax,decay,sum,
     &	        dqdall,fac,drho,tl,beta1,vnucl,rhovh,rhomu,rhoeps,
     &          eval,acc,zvnucl,utot,ekin,etot,virial,roo,q_sum
      real*8, allocatable :: rhoin(:),drofi(:),work(:),v0(:)
      allocate(rhoin(0:nr0))
      allocate(drofi(0:nr0))
      allocate(work(0:nr0))
      allocate(v0(0:nr0))
      niter=80
	q_sum=zcor(isort)
	acc=0.d0
      do li=1,nrel*lmb(isort)+1
	  if(irel.le.1) then
	    l=li-1
	  else
	    l=li/2
	  endif
        do ie=1,ntle(l,isort)
	    if(l==0) q_sum=q_sum+atoc(ie,l,isort)
	    if(l>0) q_sum=q_sum+atoc(ie,l,isort)/nrel
	  enddo
	enddo
      vhrmax=-2.d0*(z(isort)-q_sum)/rmax
	sqpi4=sqrt(4.d0*pi)
      tol=1.d-4
      tolrsq=1.d-8
      beta=0.15d0
	b=rmax/(exp(h(isort)*nr0)-1.d0)
	do ic=1,ncor(isort)
	  e_core(ic,1,isort)=-5.d0
	enddo
	do li=1,nrel*lmb(isort)+1
	  do ie=1,maxntle
	    eny(ie,li,isort,1)=-0.5d0
	  enddo
	enddo
c ------ initialize rho ------
      decay=1.d0+z(isort)/10.d0
      decay=dmin1(decay,5.d0)
      do ir=0,nr0
        roo=dexp(-decay*rofi(ir))
        rhoin(ir)=roo/sqpi4
        drofi(ir)=b*exp(h(isort)*ir)
      enddo
	do ir=0,nr0
	  work(ir)=rhoin(ir)*rofi(ir)**2*drofi(ir)
      enddo
      sum=sqpi4*dqdall(h(isort),work,nr0)
      fac=q_sum/sum
	rhoin=fac*rhoin
c -------- start self-consistency loop ------
      drho=100.d0
	nval=40
      do iter=1,niter
        tl=tolrsq
        if(iter.ge.2.and.drho.gt.2d0) tl=1.d-3
        beta1=beta
        call vcoulat(z(isort),h(isort),b,rofi,rhoin,nr0,vhrmax,v0)
        vnucl=v0(0)/sqpi4
	  do ir=0,nr0
	    work(ir)=rhoin(ir)*v0(ir)*rofi(ir)**2*drofi(ir)
        enddo
        rhovh=dqdall(h(isort),work,nr0)
	  do ir=0,nr0
	    work(ir)=rhoin(ir)*rofi(ir)*drofi(ir)
        enddo
        rhovh=rhovh-2.d0*z(isort)*sqpi4*dqdall(h(isort),work,nr0)
        call vexchat(rofi,rhoin,nr0,v0,rhomu,rhoeps,drofi,h(isort))
	  call atom(isort,v0,rofi,drofi,roatom,nr0,eval,nval)
	  if(iprn.ne.0) then
	    call prnatom(0,isort,1,acc,acc,nds)
	  endif
        do ir=0,nr0
          rhoin(ir)=beta*roatom(ir)+(1.d0-beta)*rhoin(ir)
        enddo
	  do ir=0,nr0
	    work(ir)=roatom(ir)*rofi(ir)**2*drofi(ir)
        enddo
        sum=sqpi4*dqdall(h(isort),work,nr0)
	  do ir=0,nr0
	    work(ir)=((roatom(ir)-rhoin(ir))*rofi(ir))**2*drofi(ir)
        enddo
        drho=sqrt(sqpi4*dqdall(h(isort),work,nr0))
        if(maswrk.and.iprn.ne.0) write(iun,340) iter,sum,drho
      enddo  !! over iter
#ifdef MPI
      if (nproc/=1) then
        call dgop(ecor,1,'  +',MPI_COMM_WORLD)
      endif
#endif

	zvnucl=-z(isort)*vnucl
      utot=.5d0*(rhovh+zvnucl)
      ekin=eval+ecor-rhovh-rhomu
      etot=ekin+utot+rhoeps
	virial=(utot+rhoeps)/ekin
      if(maswrk) write(iun,139) eval,ecor,vnucl,rhovh,zvnucl,utot,rhomu,
     &                          rhoeps,ekin,etot,virial
  139 format(/' sumev=',f12.5,'    sumec =',f12.5,'    vnucl=',f12.5
     .       /' rhovh=',f12.5,'    zvnucl=',f12.5,'    utot =',f12.5
     .       /' rhomu=',f12.5,'    rhoeps=',f12.5,'    ekin =',f12.5
     .       /' etot =',f12.5,'   virial =',f12.5)
     	if(nspin.eq.2) then
	  ndim=maxntle*(nrel*maxb+1)*nsort
	  call dcopy(ndim,eny,1,eny(1,1,1,2),1)
	endif
	deallocate(rhoin,drofi,work,v0)
  340 format(i5,f12.6,1p,d12.3,0p,f14.4,d14.4,f14.4,f7.2)
      end
