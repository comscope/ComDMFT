      subroutine memory_solid(maxbv,maxbj,maxw2)
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer :: m21,maxbv,maxbj,maxw2,nlmb
	nlmb=nrel*maxb+1
      m21=max0(maxb,maxpb)
	  call memory_atom
        allocate(sym(limlpb,nsort))
	  sym=.false.
        allocate(nsym(nsort))
        allocate(lmsym(limlpb,nsort))
	  lmsym=0
        if(magn.eq.2) then
          allocate(symb(limlpb,3,nsort))
	    symb=.false.
          allocate(nsymb(nsort))
          allocate(lmsymb(limlpb,3,nsort))
	    lmsymb=0
        endif
        allocate(drocor(0:maxnrad,nsort,nspin,2))
        allocate(eny(maxntle,nlmb,nsort,nspin_0))
        allocate(augm_coef(3,nlmb,nsort,nspin_0))
        allocate(fi(maxntle,nlmb,nsort,nspin))
        allocate(dfi(maxntle,nlmb,nsort,nspin))
        allocate(dny(maxntle,nlmb,nsort,nspin))
        allocate(fidot(maxntle,nlmb,nsort,nspin))
        allocate(dfidot(maxntle,nlmb,nsort,nspin))
        allocate(fidot2(maxntle,nlmb,nsort,nspin))
        allocate(dfidot2(maxntle,nlmb,nsort,nspin))
        allocate(ffsmt(2,2,maxntle,maxntle,nlmb,nsort,nspin))
        allocate(ffhmt(2,2,maxntle,maxntle,nlmb,nsort,nspin))
	  allocate(vxc0(nsort))
	  allocate(vxc1(nsort))
c	 *****************************************************
        maxbv=max(maxb,maxpb,1)
	maxbj=2*maxbv+1
      maxw2=max(maxbv,1)
	allocate(bett((maxylm+1)**2))
      call cnsbet(maxylm,bett)
      maxwig=(maxbv+1)*(2*maxbv+1)*(2*maxbv+3)/3
      maxwig_j=4*(maxbv+1)*(2*(maxbv+1)**2+1)/3
      allocate(rcostar(maxplwro/2))
      allocate(indstar(0:maxplwro/2))
      allocate(igbs(3,maxplw))
      allocate(igbsf(maxplwro))
      allocate(gbs(3,maxplw))
      allocate(ip(natom,48))
      allocate(ifia(48))
      allocate(shift(3,48))
      allocate(u(maxwig,48))
      allocate(u_opt(3,3,48))
	allocate(ag(3,3,48))
	allocate(agm(3,3,48))
      allocate(uj(maxwig_j,48))
      if(irel==2) then
        allocate(so_const(maxlfun,maxlfun,nsort))
        so_const=0.d0
      else if(irel/=2) then
        allocate(so_const_sr(maxlfun,maxlfun,nspin,nspin,natom))
        so_const_sr=0.d0
      endif
      end