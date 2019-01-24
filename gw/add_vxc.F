      subroutine add_vxc
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
      logical :: jreal,jimag
	integer :: keym
	if(ubi/='dft') return
	vxc0=0.d0
	vxc1=0.d0
	v_mt_xc=0.d0
	b_mt=0.d0
      call vexch(1,v_mt_xc,b_mt)
	v_intr_xc=(0.d0,0.d0)
	b_intr=(0.d0,0.d0)
      call vxint(1,v_intr_xc,b_intr)
      if(itest.gt.0) then
        keym=0
        if(nspin.eq.2) keym=4
	  text=' Matching : Exchange-Correlation Potential!'
        call matching(keym,v_mt_xc,v_intr_xc)
        if(nspin.eq.2) then
          keym=1
	    text=' Matching : Magnetic Field!'
          call matching(keym,v_mt_xc,v_intr_xc)
        endif
        if(magn.eq.2) then
          jreal=.false.
          jimag=.false.
          if(inv.lt.2) jreal=.true.
          if(inv.ne.1) jimag=.true.
	    text=' Matching : Magnetic Field!'
          call matchav(b_mt,b_intr,jreal,jimag)
        endif
      endif
      call timel('******* add_vxc finished ***********')
      end