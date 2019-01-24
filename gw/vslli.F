      subroutine vslli(key,key1,vhmt,vxcmt,vhint,vxcint)
c	key=0 - only Hartree part
c	key=1 - Only XC part
c	key=2 - Hartree + XC
c	key1=1 - L=0 component is taken into account
c	key1=2 - L=0 component is not taken into account
	use atom_mod
	use etot_mod
	use manager_mod
	use solid_mod
      implicit none
	integer, intent(in) :: key,key1
	real*8, intent(in) :: vhmt(maxmt/nspin_0),vxcmt(maxmt)
	complex*16, intent(in) :: vhint(nplwro),vxcint(nplwro,nspin)
	integer :: ispin,mt
	real*8, allocatable :: v(:)
	complex*16, allocatable :: vint(:)
c ---------------------------------------------------------------------
      allocate(v(maxmt/nspin_0))
      allocate(vint(nplwro))
      if(key/=0.and.magn.eq.2) call b_me
	if(key==0) then
	  v=vhmt
	  vint=vhint
	endif
	do ispin=1,nspin
	  mt=indmt(1,1,ispin)
	  if(key==1) then
	    v=vxcmt(mt:mt+maxmt/nspin_0-1)
	    vint=vxcint(:,ispin)
	  else if(key==2) then
	    v=vhmt+vxcmt(mt:mt+maxmt/nspin_0-1)
	    vint=vhint+vxcint(:,ispin)
	  endif
	  call v_me(ispin,key1,key,v)
        call warped(key,ispin,vint)
	enddo
c ---- Initial magnetic splitting ---------------------------------	  
	call split_magn(key)
c --- External magnetic field contribution ------------------------
      call h_magn_external_mt(key)	  
c ----------------------------------------------------------------	  
      call sym_me
	deallocate(v,vint)
      call timel('******* vslli finished *************')
      end