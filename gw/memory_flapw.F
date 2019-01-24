      subroutine memory_flapw
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer :: me1,nlmb
	me1=me+1
	nlmb=nrel*maxb+1
      allocate(ev_bnd(nbasmpw,nbndf,maxval(ndim3_k),nspin_0))
      allocate(z_bnd(nfun,nbndf,maxval(ndim3_k),nspin_0))
      allocate(pd_f(0:maxnrad,maxntle,nlmb,nsort,nspin))
      allocate(pd2_f(0:maxnrad,maxntle,nlmb,nsort,nspin))
      allocate(qd_f(0:maxnrad,maxntle,nlmb,nsort,nspin))
      allocate(qd2_f(0:maxnrad,maxntle,nlmb,nsort,nspin))
      allocate(vatmt(maxnsym,nsort))
	if(iter_gw+iter_qp/=0) then
        allocate(ev_qp(nbasmpw,nbndf,maxval(ndim3_k),nspin_0))
	endif
      allocate(rointr(nplwro,nspin_0))
      allocate(v_intr_h(nplwro))
      allocate(v_intr_xc(nplwro,nspin_0))
      allocate(war(nplwdif,nspin_0))
	allocate(wari(nplwdif,nspin_0))
c	if(magn.eq.2) then
	  allocate(spintr(3,nplwro))
c        if(inv.lt.2) then
	    allocate(warb(3,nplwdif))
	    warb=0.d0
c	  endif
c        if(inv.ne.1) then
	    allocate(warbi(3,nplwdif))
	    warbi=0.d0
c	  endif
	  allocate(b_intr(3,nplwro))
	  b_intr=(0.d0,0.d0)
	  allocate(b_mt(maxmtb))
c      endif
	allocate(v_ii(nplw_gw,nplw_gw,maxval(ndim3_k)))
      end
