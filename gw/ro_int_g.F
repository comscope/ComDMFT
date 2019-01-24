      SUBROUTINE ro_int_g(correlation)
	use manager_mod
	use parallel_mod
	use solid_mod
      IMPLICIT none
      logical, intent(in) :: correlation
	integer :: ispin,k,ind_k
	logical :: jreal,jimag
      DO ISPIN=1,NSPIN
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    call sumi_gw(correlation,ind_k,k,ispin)
        enddo 
        if(nproc_k/=1) then
	    call DGOP(rointr(1,ispin),2*nplwro,'  +',comm_pnt)
          if(magn.eq.2) then
	      call DGOP(spintr,6*nplwro,'  +',comm_pnt)
	    endif
	  endif
c ------- Symmetrization ----------------------------------------------
        call symscal(rointr(1,ispin))
      enddo   !!! over ispin
      if(magn.eq.2) then
        jreal=.false.
        jimag=.false.
        if(inv.lt.2) jreal=.true.
        if(inv.ne.1) jimag=.true.
        call symvec(spintr,jreal,jimag)
      endif  !!! for magn = 2 only
      END
