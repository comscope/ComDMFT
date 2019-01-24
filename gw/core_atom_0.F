      subroutine core_atom_0(v0,tol,ispin)
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
	integer, intent(in) :: ispin
	real*8, intent(in) :: tol,v0(0:nrmax)
	integer :: icore,l,li,i,lp1,it,ir,konf,nodes,kc,node,
     &	         nitec,nre,ns,ind_icor
      real*8 :: c2,sqpi4,e1,e2,val,slo,rmax,b,dff,hdif,anl,deg,ecore,
     &          convc,dl,fi2c,dqdall,de,fi00,dfi00,acc
      integer, allocatable :: ndim_ns(:),n_mpi_ns(:),ind_ns(:,:)
     	real*8, allocatable :: fc_ex(:),fc_exq(:),work(:),rocor(:),ec0(:),
     &	                     pcor0(:,:),qcor0(:,:)
      allocate(pcor0(0:nrmax,ncor(1)))
      allocate(qcor0(0:nrmax,ncor(1)))
	allocate(rocor(0:nrmax))
	allocate(fc_ex(0:nrmax),fc_exq(0:nrmax),work(0:nrmax))
	ns=ncor(1)
	allocate(ndim_ns(nproc))
	allocate(n_mpi_ns(nproc))
	call size_shift_par(ns,nproc,ndim_ns,n_mpi_ns)
	allocate(ind_ns(2,ns))
	allocate(ec0(ns))
c ---- Indexes for the parallelization --------------------------	
	icore=0
      do li=1,nrel*lmb(1)+1
	  if(irel.le.1) then
	    l=li-1
	  else
	    l=li/2
	  endif
        do konf=l+1,konfig(l,1)-1
          icore=icore+1
          ind_ns(1,icore)=konf
          ind_ns(2,icore)=li
        enddo
      enddo
	pcor0=0.d0
	qcor0=0.d0
	do ind_icor=1,ndim_ns(me+1)
	  icore=n_mpi_ns(me+1)+ind_icor
	  do ir=0,nrmax
	    pcor0(ir,icore)=pcor_atom(ir,icore,ispin)
	    qcor0(ir,icore)=qcor_atom(ir,icore,ispin)
	  enddo
	enddo
c --------------------------------------------------------------          
	fc_ex=0.d0
	fc_exq=0.d0
	c2=clight**2
	sqpi4=sqrt(4.d0*pi)
      e1=-4.d0*z(1)**2-10.d0
      e2=10.d0
      val=1.d-30
      slo=-val
      val=0.d0
      slo=-1.d-10
* define mesh parameters
      rmax=r_atom(nrmax,1)
      b=rmax/(dexp(h(1)*nrmax-h(1))-1.d0)
      hdif=1.d0/h(1)
c ------ start loop over core states --------------
	icore=0
	rocor=0.d0
	ec0=0.d0
	cor_norm(:,ispin,1)=0.d0
	r_nre_core(:,ispin,1)=0.d0
	psi_nre(:,ispin,1)=0.d0
	do ind_icor=1,ndim_ns(me+1)
	  icore=n_mpi_ns(me+1)+ind_icor
        konf=ind_ns(1,icore)
        li=ind_ns(2,icore)
	  if(irel.le.1) then
	    l=li-1
	    anl=dfloat(l+l+1)
	  else
	    l=li/2
	    anl=dfloat(((li+1)/2)*2)
	    i=li-2*l
	    if(i.eq.0) i=-1
	  endif
        lp1=l+1
        deg=2.d0*anl/nspin/nrel
        nodes=konf-lp1
        ecore=e_core(icore,ispin,1)
	  it=0
	  dff=1.d0
	  nre=nrmax
	  do while(it<10.and.dff>1.d-10)
	    call f_ex_atom(fc_ex,fc_exq,pcor0(0,icore),qcor0(0,icore),
     &	                 ispin,li)
	    call rad_eq(e1,e2,ecore,tol,z(1),l,i,nodes,val,slo,v0(0),fi2c,
     &	              r_atom(0,1),nrmax,pcor0(0,icore),qcor0(0,icore),
     &                fc_ex,fi00,dfi00,kc,nitec,convc,nre,fc_exq,
     &                dr_atom,1,h(1),dl,work,acc,node,irel)
	    it=it+1
          dff=abs(e_core(icore,ispin,1)-ecore)
          e_core(icore,ispin,1)=ecore
          cor_norm(icore,ispin,1)=fi2c
          r_nre_core(icore,ispin,1)=r_atom(nre,1)
          psi_nre(icore,ispin,1)=pcor0(nre,icore)/r_atom(nre,1)
          ec0(icore)=ecore
        enddo
	  ecor=ecor+ecore*deg
	  do ir=1,nrmax
	    rocor(ir)=rocor(ir)+deg*(pcor0(ir,icore)/r_atom(ir,1))**2
	  enddo
	  if(irel.ge.1) then
	    do ir=1,nrmax
	      rocor(ir)=rocor(ir)+deg*(qcor0(ir,icore)/r_atom(ir,1))**2/c2
	    enddo
	  endif
	  work(0)=0.d0
	  do ir=1,nrmax
	    work(ir)=pcor0(ir,icore)*fc_ex(ir)*dr_atom(ir)
	  enddo
	  if(irel.ge.1) then
	    do ir=1,nrmax
	      work(ir)=work(ir)+qcor0(ir,icore)*fc_exq(ir)*dr_atom(ir)/c2
	    enddo
	  endif
        gx_sx_cor=gx_sx_cor+dqdall(h(1),work,nrmax)*deg
	  de=ecore-chem_pot
      enddo
c -------------------------------------------------------------------
      if(nproc/=1.and.ns/=0) then
#ifdef MPI
	  call dgop(rocor,nrmax+1,'  +',MPI_COMM_WORLD)
	  call dgop(ec0,ns,'  +',MPI_COMM_WORLD)
	  call dgop(pcor0,ns*(nrmax+1),'  +',MPI_COMM_WORLD)
	  call dgop(qcor0,ns*(nrmax+1),'  +',MPI_COMM_WORLD)
	  call dgop(cor_norm(1,ispin,1),ncormax,'  +',MPI_COMM_WORLD)
	  call dgop(r_nre_core(1,ispin,1),ncormax,'  +',MPI_COMM_WORLD)
	  call dgop(psi_nre(1,ispin,1),ncormax,'  +',MPI_COMM_WORLD)
#endif	  
      endif
      e_core(1:ns,ispin,1)=ec0
	do icore=1,ncor(1)
	  do ir=0,nrmax
	    pcor_atom(ir,icore,ispin)=pcor0(ir,icore)
	    qcor_atom(ir,icore,ispin)=qcor0(ir,icore)
	  enddo
	enddo
c -------------------------------------------------------------------      
	work(0)=0.d0
      do ir=1,nrmax
        work(ir)=rocor(ir)*(v0(ir)-2.d0*z(1)/r_atom(ir,1))
     &	                               *dr_atom(ir)*r_atom(ir,1)**2
      enddo
      gx_vh_cor=gx_vh_cor+dqdall(h(1),work,nrmax)
      do ir=1,nrmax
        ro_core(ir,1,ispin)=rocor(ir)/sqpi4
      enddo
	deallocate(pcor0,qcor0)
	deallocate(fc_ex,fc_exq,work,rocor,ndim_ns,n_mpi_ns,ind_ns,ec0)
      end
      
      
      subroutine f_ex_atom(f_ex,f_exq,pp,qq,ispin,li)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ispin,li
	real*8, intent(in) :: pp(0:nrmax),qq(0:nrmax)
	real*8, intent(out) :: f_ex(0:nrmax),f_exq(0:nrmax)
	integer :: li1,lf1,l1,lf2,li2,l2,l,ic1,ir,ll2,je,ie
	real*8 :: pi8,deg,fc
	real*8, allocatable :: pval(:),qval(:),pp1(:),f(:),g(:),qq1(:),
     &	                     tmp(:),work(:),ff(:)
	allocate(tmp(0:nrmax),pp1(0:nrmax),f(0:nrmax),g(0:nrmax),
     &	       qval(0:nrmax),qq1(0:nrmax),work(0:nrmax),pval(0:nrmax))
      allocate(ff(0:nrmax))
      pi8=8.d0*pi
	f_ex(0:nrmax)=0.d0
	f_exq(0:nrmax)=0.d0
	if(irel.le.1) then
	  l=li-1
	  deg=l+l+1.d0
	else
	  l=li/2
	  deg=dfloat(((li+1)/2)*2)
	endif
	fc=pi8/deg
c ------ core contribution --------------
	do ic1=1,ncor(1)
	  li1=li_cor(ic1,1)
	  if(irel.le.1) then
	    l1=li1-1
	  else
	    l1=li1/2
	  endif
	  pp1=pcor_atom(:,ic1,ispin)
	  if(irel.ge.1) qq1=qcor_atom(:,ic1,ispin)
	  do l2=0,min(l+l1,max_rl)
	    call var_int(dr_atom(0),pp(0),qq(0),pp1(0),qq1(0),f(0),g(0),
     &	               tmp(0),r_atom(0,1),rl_atom(0,l2),work(0),
     &                 nrmax+1,h(1))
	    f=-tmp*d_lll(l2,li,li1)**2*fc
	    do ir=0,nrmax
	      f_ex(ir)=f_ex(ir)+f(ir)*pp1(ir)
	    enddo
	    if(irel.ge.1) then
	      do ir=0,nrmax
	        f_exq(ir)=f_exq(ir)+f(ir)*qq1(ir)
	      enddo
	    endif
	  enddo   !! over l2
	enddo  !! over ic1
c ------ Valence contribution -----------------------------------------
      do li2=1,nrel*lmb(1)+1
	  if(irel.le.1) then
	    ll2=li2-1
	  else
	    ll2=li2/2
	  endif
        do je=1,ntle(ll2,1)
          lf2=indfun0_atom(je,li2)
	    pval=pval_atom(:,lf2,ispin)
	    if(irel.ge.1) qval=qval_atom(:,lf2,ispin)
	    ff=0.d0
	    do l2=0,min(l+ll2,max_rl)
	      call var_int(dr_atom(0),pval(0),qval(0),pp(0),qq(0),
     &		  	         f(0),g(0),tmp(0),r_atom(0,1),rl_atom(0,l2),
     &                   work(0),nrmax+1,h(1))
	      do ir=0,nrmax
	        ff(ir)=ff(ir)+tmp(ir)*d_lll(l2,li,li2)**2
	      enddo
	    enddo
          do ie=1,ntle(ll2,1)
            lf1=indfun0_atom(ie,li2)
	      f=pval_atom(:,lf1,ispin)
	      if(irel>=1) g=qval_atom(:,lf1,ispin)
	      do ir=0,nrmax
	        f_ex(ir)=f_ex(ir)+ff(ir)*f(ir)
     &	                       *g_beta_atom(ie,je,li2,ispin)*fc
	      enddo
	      if(irel>=1) then
	        do ir=0,nrmax
	          f_exq(ir)=f_exq(ir)+ff(ir)*g(ir)*fc
     &	                           *g_beta_atom(ie,je,li2,ispin)
 	        enddo
	      endif
	    enddo   !! over lf1
	  enddo   !! over l2
	enddo  !! over lf2
	deallocate(tmp,pp1,f,g,qval,qq1,work,pval,ff)
      end