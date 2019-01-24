      subroutine e_x_cor(isort,ispin,ex)
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
	integer, intent(in) :: isort,ispin
	real*8, intent(out) :: ex
	integer :: icore,l,li,nr0,i,lp1,ir,konf,nodes,ns,ind_icor,mt
      real*8 :: c2,anl,dqdall
      integer, allocatable :: ndim_ns(:),n_mpi_ns(:),ind_ns(:,:)
     	real*8, allocatable :: fc_ex(:),fc_exq(:),work(:),pcor0(:,:),
     &	                     qcor0(:,:)
     	nr0=nrad(isort)
	allocate(fc_ex(0:nr0),fc_exq(0:nr0),work(0:nr0))
	ns=ncor(isort)
	allocate(ndim_ns(nproc))
	allocate(n_mpi_ns(nproc))
	call size_shift_par(ns,nproc,ndim_ns,n_mpi_ns)
	allocate(ind_ns(2,ns))
c ---- Indexes for the parallelization --------------------------	
	icore=0
      do li=1,nrel*lmb(isort)+1
	  if(irel.le.1) then
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
	allocate(pcor0(0:nr0,ndim_ns(me+1)))
	allocate(qcor0(0:nr0,ndim_ns(me+1)))
	pcor0=0.d0
	qcor0=0.d0
	do ind_icor=1,ndim_ns(me+1)
	  icore=n_mpi_ns(me+1)+ind_icor
	  mt=indcor(icore,ispin,isort)
	  do ir=0,nrad(isort)
	    pcor0(ir,ind_icor)=pcor_old(mt+ir)
	    qcor0(ir,ind_icor)=qcor_old(mt+ir)
	  enddo
	enddo
c --------------------------------------------------------------          
	fc_ex=0.d0
	fc_exq=0.d0
	c2=clight**2
	ex=0.d0
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
        nodes=konf-lp1
	  call f_ex_new(fc_ex,fc_exq,isort,pcor0(0,ind_icor),
     &                qcor0(0,ind_icor),ispin,li,nr0)
	  work(0)=0.d0
	  do ir=1,nr0
	    work(ir)=pcor0(ir,ind_icor)*fc_ex(ir)*dr(ir,isort)
	  enddo
	  if(irel.ge.1) then
	    do ir=1,nr0
	      work(ir)=work(ir)+qcor0(ir,ind_icor)*fc_exq(ir)*dr(ir,isort)
     &	                                        /c2
	    enddo
	  endif
        ex=ex+dqdall(h(isort),work,nr0)*anl
      enddo
c -------------------------------------------------------------------
#ifdef MPI
      if(nproc/=1.and.ns/=0) call dgop(ex,1,'  +',MPI_COMM_WORLD)
#endif      
	deallocate(fc_ex,fc_exq,work,ndim_ns,n_mpi_ns,ind_ns,pcor0,qcor0)
      end