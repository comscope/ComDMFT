      subroutine tos_lda(x,ttos,tdos)
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer :: ispin,k,ibnd,n
      real*8 :: const,st,sd,ffd,fermi_dirac,x
      real*8 :: ttos ! sum of Fermi-Dirac occupation numbers
                     ! should equal total number of electrons
      real*8 :: tdos ! derivative of ttos wrt Fermi energy
	const=2.d0/dfloat(nspin*nrel)
      ttos=0.d0
      tdos=0.d0
      do ispin=1,nspin
        do k=1,npnt
          n=n_bnd(k,ispin)
          st=0.d0
          sd=0.d0
          do ibnd=1,n
            ffd=fermi_dirac(e_bnd(ibnd,k,ispin)-x)
            st=st+ffd
            sd=sd+betta_t*ffd*(1.d0-ffd)
          enddo
          ttos=ttos+wgt(k)*st*const
          tdos=tdos+wgt(k)*sd*const
        enddo  !!! over k
      enddo  !!! over ispin
      if(maswrk) write(iun,'(3(1x,f13.8))')x,ttos,tdos
      end
