      subroutine rofull_atom
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use units_mod
	use solid_mod
      implicit none
	integer :: li,ir,ival,l,ie,ispin,ie1,ival1
      real*8 :: c2,sqpi4,deg,de,fc,deg1
      ro_atom=0.d0
      do ispin=1,nspin
        ro_atom(:,1)=ro_atom(:,1)+ro_core(:,1,ispin)
      enddo
	c2=clight**2
	sqpi4=sqrt(4.d0*pi)
c ------ loop over valence states -------
      fc=-2.d0/sqpi4/nrel/nspin
      do ispin=1,nspin
        do li=1,nrel*lmb(1)+1
	    if(irel.le.1) then
	      l=li-1
	      deg=l+l+1.d0
	    else
	      l=li/2
	      deg=dfloat(((li+1)/2)*2)
	    endif
          do ie1=1,ntle(l,1)
	      ival1=indfun0_atom(ie1,li)
            do ie=1,ntle(l,1)
	        ival=indfun0_atom(ie,li)
	        deg1=fc*deg*g_beta_atom(ie,ie1,li,ispin)
	        do ir=1,nrmax
	          de=pval_atom(ir,ival,ispin)*pval_atom(ir,ival1,ispin)
	          if(irel>=1) de=de+qval_atom(ir,ival,ispin)
     &	                         *qval_atom(ir,ival1,ispin)/c2
	          ro_atom(ir,1)=ro_atom(ir,1)+deg1*de/r_atom(ir,1)**2
	        enddo
	      enddo
	    enddo
	  enddo
	enddo
	ro_atom(:,1)=admix*ro_atom(:,1)+(1.d0-admix)*ro_atom_old
	ro_atom_old=ro_atom(:,1)
      end