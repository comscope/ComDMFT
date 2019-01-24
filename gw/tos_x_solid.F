      subroutine tos_x_solid(x,ttos,tdos)
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	real*8, intent(in) :: x
	real*8, intent(inout) :: ttos,tdos
	integer :: me1,ispin,k,i,n
	real*8 :: const,st,sd,de,ffd,dffd,fermi_dirac
	me1=me+1
	const=2.d0/dfloat(nspin*nrel)
      do ispin=1,nspin
        do k=1,npnt
	    n=n_bnd(k,ispin)
          st=0.d0
          sd=0.d0
	    do i=1,n
            de=e_bnd(i,k,ispin)-x
            ffd=fermi_dirac(de)
            dffd=betta_t*ffd*(1.d0-ffd)
            st=st+ffd
            sd=sd+dffd
          enddo
          ttos=ttos+wgt(k)*st*const
          tdos=tdos+wgt(k)*sd*const
        enddo  !!! over k
      enddo  !!! over ispin
      if(maswrk.and.iprn.ge.1) write(iun,'(3(1x,f13.8))')x,ttos,tdos
      end