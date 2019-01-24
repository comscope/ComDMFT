      real*8 function sum_mats_sigx_g_loc(n,sigx,g)
      use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n
	real*8, intent(in) :: sigx(maxel_loc,maxel_loc,nspin),
     &	                    g(maxel_loc,maxel_loc,2,ndim3_tau,nspin)
	integer :: j,i,ispin
	real*8 :: e
	e=0.d0
      if(me3_tau==0) then
        do ispin=1,nspin
	    do j=1,n
	      do i=1,n
	        e=e+sigx(i,j,ispin)*g(j,i,2,1,ispin)
	      enddo
	    enddo
        enddo
        e=-0.5d0*e
      endif
	if(nproc_tau/=1) call dgop(e,1,'  +',comm_pnk)
	sum_mats_sigx_g_loc=e*2.d0/nspin
      end
      
      
      real*8 function sum_mats_sigc_g_loc(n,sig,g)
      use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n
	complex*16, intent(in) :: sig(maxel_loc,maxel_loc,ndim3_omega,
     &	                            nspin),
     &	                        g(maxel_loc,maxel_loc,ndim3_omega,
     &                            nspin)
	integer :: i_omega,j,i,ind_omega,ispin
	real*8 :: de,e
	e=0.d0
      do ispin=1,nspin
	  do ind_omega=1,ndim3_omega
	    i_omega=me3_tau*ndim3_omega+ind_omega-1
	    de=0.d0
	    do j=1,n
	      do i=1,n
	        de=de+sig(i,j,ind_omega,ispin)*g(j,i,ind_omega,ispin)
	      enddo
	    enddo
          e=e+et_sum_sc_gc(i_omega)*de
        enddo
      enddo
	if(nproc_tau/=1) call dgop(e,1,'  +',comm_pnk)
	sum_mats_sigc_g_loc=e*2.d0/nspin
      end