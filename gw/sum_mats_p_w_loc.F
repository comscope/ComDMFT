      subroutine sum_mats_p_w_loc(n,np,p_tau,p_nu,upb,ws_nu,sigx,ex,
     &                            ec)
      use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n,np
	real*8, intent(in) :: p_tau(np,np,ndim3_tau),p_nu(np,np,ndim3_nu),
     &	                    ws_nu(np,np,ndim3_nu),upb(np,np),
     &                      sigx(n,n,nspin)
      real*8, intent(out) :: ex,ec
	integer :: i_nu,j,i,ind_nu,ispin
	real*8 :: de
	ex=0.d0
      if(me3_tau==0) then
        do ispin=1,nspin
          do i=1,n
            ex=ex+sigx(i,i,ispin)
          enddo
        enddo
        ex=2.d0*ex/nspin
	  do j=1,np
	    do i=1,np
	      ex=ex-p_tau(i,j,1)*upb(j,i)
	    enddo
        enddo
        ex=0.5d0*ex
      endif
	if(nproc_tau/=1) call dgop(ex,1,'  +',comm_pnk)
	ec=0.d0
	do ind_nu=1,ndim3_nu
	  i_nu=me3_tau*ndim3_nu+ind_nu-1
	  de=0.d0
	  do j=1,np
	    do i=1,np
	      de=de-p_nu(i,j,ind_nu)*ws_nu(j,i,ind_nu)
	    enddo
	  enddo
        ec=ec+et_sum_nu(i_nu)*de
      enddo
	if(nproc_tau/=1) call dgop(ec,1,'  +',comm_pnk)
      end