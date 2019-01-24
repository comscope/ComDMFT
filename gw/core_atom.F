      subroutine core_atom
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
	integer :: ispin,li,icor,ir,ii,jcor,lj,je,i
	real*8 :: sqpi4,c2,ort(6),dqdall,nrm
	real*8, allocatable :: work(:)
	if(ncor_tot==0) return
	c2=clight**2
	sqpi4=sqrt(4.d0*pi)
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c      do i=1,9
	ecor=0.d0
	gx_vh_cor=0.d0
	gx_sx_cor=0.d0
	gx_vxc_cor=0.d0
	vhr=vhr/sqpi4
	do ispin=1,nspin
		call core_atom_0(vhr(0),1.d-14,ispin)
	enddo
	vhr=vhr*sqpi4
c	enddo
	if(nproc/=1) then
#ifdef MPI
	  call dgop(ecor,1,'  +',MPI_COMM_WORLD)
	  call dgop(gx_sx_cor,1,'  +',MPI_COMM_WORLD)
#endif	  
      endif
      allocate(work(0:nrmax))
      if(maswrk) then
	  do ispin=1,nspin
	    if(ncor(1).ne.0) then
	      write(iun,*)' Spin = ',ispin
		    write(iun,*)' Core states'
	      write(iun,201)
          endif
	    do icor=1,ncor(1)
	      li=li_cor(icor,1)
c ----- MT-check check ---------------------------
	      work(0)=0.d0
	      do ir=0,nrmax
	        work(ir)=pcor_atom(ir,icor,ispin)**2*dr_atom(ir)
	      enddo
	      if(irel.ge.1) then
	        do ir=0,nrmax
	          work(ir)=work(ir)+qcor_atom(ir,icor,ispin)**2
     &	                         *dr_atom(ir)/c2
	        enddo
	      endif
            nrm=dqdall(h(1),work,nrmax)
	      write(iun,202)li,e_core(icor,ispin,1),
     &	                  pcor_atom(nrmax,icor,ispin)/r_atom(nrmax,1),
     &	                  psi_nre(icor,ispin,1),
     &	                  r_nre_core(icor,ispin,1),
     &                    cor_norm(icor,ispin,1)-1.d0,nrm
          enddo
		  write(iun,*)' Orthogonality check'
	    do icor=1,ncor(1)
	      li=li_cor(icor,1)
c ----- Orthogonality check ---------------------------
            ii=0
	      do jcor=icor+1,ncor(1)
	        lj=li_cor(jcor,1)
	        if(lj==li) then
	          ii=ii+1
	          je=indcor(jcor,ispin,1)
	          work(0)=0.d0
	          do ir=0,nrmax
	            work(ir)=pcor_atom(ir,icor,ispin)
     &	                  *pcor_atom(ir,jcor,ispin)*dr_atom(ir)
	          enddo
	          if(irel.ge.1) then
	            do ir=0,nrmax
	              work(ir)=work(ir)+qcor_atom(ir,icor,ispin)
     &	                             *qcor_atom(ir,jcor,ispin)
     &	                             *dr_atom(ir)/c2
	            enddo
	          endif
                ort(ii)=dqdall(h(1),work,nrmax)
              endif
            enddo
c -------------------------------------------------- 
	      write(iun,203)li,(ort(i),i=1,ii)
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