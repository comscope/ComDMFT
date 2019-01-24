      subroutine cor_new(v0,tol,rofi,drofi,isort,nr0,pcor0,qcor0,isp,vh,
     &                   vx,key1)
c	key1=1 - Exchange is included
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
	integer, intent(in) :: key1,nr0,isp,isort
	real*8, intent(in) :: v0(0:nr0),tol,rofi(0:nr0),drofi(0:nr0),
     &	                    vh(0:nr0),vx(0:nr0)
	real*8, intent(inout) :: pcor0(0:nrmax,ncor(isort)),
     &	                       qcor0(0:nrmax,ncor(isort))
	integer :: icore,l,li,i,lp1,it,ir,konf,nodes,kc,nitec,nre,ns,
     &	         ind_icor,mt,node
      real*8 :: c2,sqpi4,e1,e2,val,slo,rmax,dff,hdif,anl,deg,ecore,
     &          convc,dl,fi2c,dqdall,de,fi00,dfi00,acc
      integer, allocatable :: ndim_ns(:),n_mpi_ns(:),ind_ns(:,:)
     	real*8, allocatable :: fc_ex(:),fc_exq(:),work(:),rocor(:),ec0(:)
	allocate(rocor(0:nr0))
	allocate(fc_ex(0:nr0),fc_exq(0:nr0),work(0:nr0))
	ns=ncor(isort)
	allocate(ndim_ns(nproc))
	allocate(n_mpi_ns(nproc))
	call size_shift_par(ns,nproc,ndim_ns,n_mpi_ns)
	allocate(ind_ns(2,ns))
	allocate(ec0(ns))
c ---- Indexes for the parallelization --------------------------	
	icore=0
      do li=1,nrel_core*lmb(isort)+1
	  if(irel_core.le.1) then
	    l=li-1
	  else
	    l=li/2
	  endif
        do konf=l+1,konfig(l,isort)-1
          icore=icore+1
          ind_ns(1,icore)=konf
          ind_ns(2,icore)=li
        enddo
      enddo
	pcor0=0.d0
	qcor0=0.d0
	do ind_icor=1,ndim_ns(me+1)
	  icore=n_mpi_ns(me+1)+ind_icor
	  mt=indcor(icore,isp,isort)
	  do ir=0,nrad(isort)
	    pcor0(ir,icore)=pcor(mt+ir)
	    qcor0(ir,icore)=qcor(mt+ir)
	  enddo
	enddo
c --------------------------------------------------------------          
	fc_ex=0.d0
	fc_exq=0.d0
	c2=clight**2
	sqpi4=sqrt(4.d0*pi)
      e1=-4.d0*z(isort)**2-10.d0
      e2=10.d0
      if(core_mt) val=0.d0
      if(.not.core_mt) val=1.d-10
      slo=-1.d-10
* define mesh parameters
      rmax=rofi(nr0)
      hdif=1.d0/h(isort)
c ------ start loop over core states --------------
	icore=0
	rocor=0.d0
	ec0=0.d0
	cor_norm(:,isp,isort)=0.d0
	r_nre_core(:,isp,isort)=0.d0
	psi_nre(:,isp,isort)=0.d0
	do ind_icor=1,ndim_ns(me+1)
	  icore=n_mpi_ns(me+1)+ind_icor
        konf=ind_ns(1,icore)
        li=ind_ns(2,icore)
	  if(irel_core.le.1) then
	    l=li-1
	    anl=dfloat(l+l+1)
	  else
	    l=li/2
	    anl=dfloat(((li+1)/2)*2)
	    i=li-2*l
	    if(i.eq.0) i=-1
	  endif
        lp1=l+1
        deg=2.d0*anl/nspin/nrel_core
        nodes=konf-lp1
        ecore=e_core(icore,isp,isort)
	  it=0
	  dff=1.d0
	  nre=nr0
	  do while(it<10) !.and.dff>1.d-10)
	    if(key1==1) call f_ex_new(fc_ex,fc_exq,isort,pcor0(0,icore),
     &                              qcor0(0,icore),isp,li,nr0)
	    call rad_eq(e1,e2,ecore,tol,z(isort),l,i,nodes,val,slo,
     &                v0(0),fi2c,rofi,nr0,pcor0(0,icore),
     &                qcor0(0,icore),fc_ex,fi00,dfi00,kc,nitec,convc,
     &                nre,fc_exq,drofi,key1,h(isort),dl,work,
     &                acc,node,irel_core)
	    it=it+1
          dff=abs(e_core(icore,isp,isort)-ecore)
          e_core(icore,isp,isort)=ecore
          if(ubi/='dft') cor_norm(icore,isp,isort)=fi2c
          if(ubi=='dft') cor_norm(icore,isp,isort)=1.d0
          r_nre_core(icore,isp,isort)=rofi(nre)
          psi_nre(icore,isp,isort)=pcor0(nre,icore)/rofi(nre)
          ec0(icore)=ecore
        enddo
	  ecor=ecor+ecore*deg*nhsort(isort)
	  do ir=1,nr0
	    rocor(ir)=rocor(ir)+deg*(pcor0(ir,icore)/rofi(ir))**2
	  enddo
	  if(irel_core.ge.1) then
	    do ir=1,nr0
	      rocor(ir)=rocor(ir)+deg*(qcor0(ir,icore)/rofi(ir))**2/c2
	    enddo
	  endif
	  if(key1==1) then
	    work(0)=0.d0
	    do ir=1,nr0
	      work(ir)=pcor0(ir,icore)*fc_ex(ir)*drofi(ir)
	    enddo
	    if(irel_core.ge.1) then
	      do ir=1,nr0
	        work(ir)=work(ir)+qcor0(ir,icore)*fc_exq(ir)*drofi(ir)/c2
	      enddo
	    endif
          gx_sx_cor=gx_sx_cor+dqdall(h(isort),work,nr0)*deg
     &                                                 *nhsort(isort)
	  endif
	  de=ecore-chem_pot
      enddo
c -------------------------------------------------------------------
      if(nproc/=1.and.ns/=0) then
#ifdef MPI
	  call dgop(rocor,nr0+1,'  +',MPI_COMM_WORLD)
	  call dgop(ec0,ns,'  +',MPI_COMM_WORLD)
	  call dgop(pcor0,ns*(nrmax+1),'  +',MPI_COMM_WORLD)
	  call dgop(qcor0,ns*(nrmax+1),'  +',MPI_COMM_WORLD)
	  call dgop(cor_norm(1,isp,isort),ncormax,'  +',MPI_COMM_WORLD)
	  call dgop(r_nre_core(1,isp,isort),ncormax,'  +',MPI_COMM_WORLD)
	  call dgop(psi_nre(1,isp,isort),ncormax,'  +',MPI_COMM_WORLD)
#endif	  
      endif
      e_core(1:ns,isp,isort)=ec0
	do icore=1,ncor(isort)
	  mt=indcor(icore,isp,isort)
	  do ir=0,nrad(isort)
	    pcor(mt+ir)=pcor0(ir,icore)
	    qcor(mt+ir)=qcor0(ir,icore)
	  enddo
	enddo
c -------------------------------------------------------------------      
	work(0)=0.d0
      do ir=1,nr0
        work(ir)=rocor(ir)*(vh(ir)/sqpi4-2.d0*z(isort)/rofi(ir))
     &	                                      *drofi(ir)*rofi(ir)**2
      enddo
      gx_vh_cor=gx_vh_cor+dqdall(h(isort),work,nr0)*nhsort(isort)
	if(key1==0) then 
	  work(0)=0.d0
        do ir=1,nr0
          work(ir)=rocor(ir)*vx(ir)/sqpi4*drofi(ir)*rofi(ir)**2
        enddo
        gx_vxc_cor=gx_vxc_cor+dqdall(h(isort),work,nr0)*nhsort(isort)
	endif
      do ir=1,nr0
        ro_core(ir,isort,isp)=rocor(ir)/sqpi4
      enddo
	deallocate(fc_ex,fc_exq,work,rocor,ndim_ns,n_mpi_ns,ind_ns,ec0)
      end