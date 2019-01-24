      subroutine core_all
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
	integer :: isort,ispin,mt,mt0,key1,li,ie,icor,ir,nr0,ii,jcor,
     &	         lj,je,i,mt1
	real*8 :: sqpi4,conv,betta,c2,ort(6),dqdall,nrm,dd,a,b
	real*8, allocatable :: v0(:),vh(:),vx(:),rc(:),drc(:),pcor0(:,:),
     &	                     qcor0(:,:),work(:)
	if(ncor_tot==0) return
	c2=clight**2
	sqpi4=sqrt(4.d0*pi)
	pcor_old=pcor
	qcor_old=qcor
	allocate(t_x(nrel_core*maxl_cor+1,nrel_core*maxl_cor+1,
     &	           0:2*maxl_cor),
     &         t1_x(0:maxb+maxb,maxlfun,maxlfun,nrel_core*maxl_cor+1))
      allocate(rc(0:nrmax),drc(0:nrmax))
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
      key1=0
      if(ubi/='dft') key1=1
	allocate(v0(0:nrmax))
	allocate(vh(0:nrmax))
	allocate(vx(0:nrmax))
	ecor=0.d0
	gx_vh_cor=0.d0
	gx_sx_cor=0.d0
	gx_vxc_cor=0.d0
	do isort=1,nsort
	  if(ncor(isort)==0) cycle
        mt0=indmt(1,isort,1)
        mt1=indmt(1,isort,nspin)
        betta=smt(isort)/(exp(h(isort)*nrad(isort))-1.d0)
        allocate(pcor0(0:nrmax,ncor(isort)))
        allocate(qcor0(0:nrmax,ncor(isort)))
	  do ispin=1,nspin
          mt=indmt(1,isort,ispin)
	    if(key1==1) call t_t1_x(isort,ispin,0)
	    vh(0:nrad(isort))=v_mt_h(mt0:mt0+nrad(isort))
	    v0(0:nrad(isort))=vh(0:nrad(isort))
	    if(key1==0) then
	      vx(0:nrad(isort))=v_mt_xc(mt:mt+nrad(isort))
	      v0(0:nrad(isort))=v0(0:nrad(isort))+vx(0:nrad(isort))
	    endif
	    v0(0:nrad(isort))=v0(0:nrad(isort))/sqpi4
	    rc(0:nrad(isort))=r(0:nrad(isort),isort)
	    drc(0:nrad(isort))=dr(0:nrad(isort),isort)
	    a=vh(nrad(isort))-2*z(isort)/rc(nrad(isort))*sqpi4
	    b=vh(nrad(isort)-1)-2*z(isort)/rc(nrad(isort)-1)*sqpi4
	    dd=(a-b)/(rc(nrad(isort))-rc(nrad(isort)-1))
          if(dd<0.d0) dd=0.d0
	    do ir=nrad(isort)+1,nr(isort)
            rc(ir)=betta*(exp(h(isort)*ir)-1.d0)
            drc(ir)=betta*exp(h(isort)*ir)
c	      vh(ir)=vh(nrad(isort))+dd*(rc(ir)-rc(nrad(isort)))
	      vh(ir)=a+2*z(isort)/rc(ir)*sqpi4
	      v0(ir)=vh(ir)/sqpi4
	      if(key1==0) then
	        vx(ir)=0.d0  !vx(nrad(isort))
	        v0(ir)=(vh(ir)+vx(ir))/sqpi4
	      endif
	    enddo
	    if(core_mt) nr0=nrad(isort)
	    if(.not.core_mt) nr0=nr(isort)
		  call cor_new(v0(0),1.d-14,rc,drc,isort,nr0,pcor0,qcor0,
     &		           ispin,vh,vx,key1)
	  enddo
	  deallocate(pcor0,qcor0)
	enddo
	if(nproc/=1) then
#ifdef MPI
	  call dgop(ecor,1,'  +',MPI_COMM_WORLD)
	  call dgop(gx_sx_cor,1,'  +',MPI_COMM_WORLD)
#endif	  
      endif
	conv=maxval(abs(pcor_old-pcor))
	if(maswrk) write(iun,*)' Core conv. ',conv
	deallocate(v0,vh,vx,t_x,t1_x,rc,drc)
      allocate(work(0:nrmax))
      if(maswrk) then
	  do isort=1,nsort
	    do ispin=1,nspin
	      if(ncor(isort).ne.0) then
	        write(iun,*)' Spin = ',ispin
		      write(iun,*)' Core states'
	        write(iun,201)
            endif
	      do icor=1,ncor(isort)
	        li=li_cor(icor,isort)
	        ie=indcor(icor,ispin,isort)
c ----- MT-check check ---------------------------
	        work(0)=0.d0
	        do ir=0,nrad(isort)
	          work(ir)=pcor(ie+ir)*pcor(ie+ir)*dr(ir,isort)
	        enddo
	        if(irel_core.ge.1) then
	          do ir=0,nrad(isort)
	            work(ir)=work(ir)+qcor(ie+ir)*qcor(ie+ir)
     &	                           *dr(ir,isort)/c2
	          enddo
	        endif
              nrm=dqdall(h(isort),work,nrad(isort))
	        write(iun,202)li,e_core(icor,ispin,isort),
     &	                    pcor(ie+nrad(isort))/smt(isort),
     &	                    psi_nre(icor,ispin,isort),
     &	                    r_nre_core(icor,ispin,isort),
     &                      cor_norm(icor,ispin,isort)-1.d0,nrm
            enddo
		    write(iun,*)' Orthogonality check'
	      do icor=1,ncor(isort)
	        li=li_cor(icor,isort)
	        ie=indcor(icor,ispin,isort)
c ----- Orthogonality check ---------------------------
              ii=0
	        do jcor=icor+1,ncor(isort)
	          lj=li_cor(jcor,isort)
	          if(lj==li) then
	            ii=ii+1
	            je=indcor(jcor,ispin,isort)
	            work(0)=0.d0
	            do ir=0,nrad(isort)
	              work(ir)=pcor(ie+ir)*pcor(je+ir)*dr(ir,isort)
	            enddo
	            if(irel_core.ge.1) then
	              do ir=0,nrad(isort)
	                work(ir)=work(ir)+qcor(ie+ir)*qcor(je+ir)
     &	                               *dr(ir,isort)/c2
	              enddo
	            endif
                  ort(ii)=dqdall(h(isort),work,nrad(isort))
                endif
              enddo
c -------------------------------------------------- 
	        write(iun,203)li,(ort(i),i=1,ii)
            enddo
          enddo
	  enddo
        call timel('** CORE finished *******************')
      endif
      deallocate(work)
201	format(' li',6x,'E',12x,'Psi_MT',5x,' Psi_NRE ',2x,' R_NRE ',2x,
     &       'NORM-1',6x,'Norm_MT')
202	format(i3,g14.7,2(1x,e12.5),f7.3,e12.5,1x,f9.6)
203	format(i3,6e12.5)
      end