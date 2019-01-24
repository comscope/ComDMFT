      subroutine memory_atom
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer :: m21,lllmax,i1,i2,nfac,ikey,icg
	integer, allocatable :: ifcg(:),ijcg(:)
      m21=max0(maxb,maxpb)
	lllmax=23
      if(m21.gt.lllmax) then
        if(maswrk) write(iun,*) 'm21 > lllmax'
        call ending
      endif
	allocate(ifcg(0:lllmax),ijcg(0:lllmax))
	ifcg=(/2,16,97,393,1188,2958,6409,12537,22662,38492,62153,
     &96257,143928,208874,295409,408529,553930,738088,968273,1252633,
     &1600204,2020998,2526009,3127305/)
      ijcg=(/2,11,46,137,326,667,1226,2081,3322,5051,7382,
     &10441,14366,19307,25426,32897,41906,52651,65342,80201,97462,
     &117371,140186,166177/)
      i1=ifcg(m21)
      i2=ijcg(m21)
	deallocate(ifcg,ijcg)
      allocate(cg(i1))
      allocate(indxcg(i2))
      allocate(jcg(i1))
      maxylm=max(2*maxb,maxpb,20)
	max_rl=2*maxb
      nfac=100
      call scg(m21,cg,indxcg,jcg)
      if(irel.eq.2.or.irel_core==2) then
        call SGAUNT
        call sckg2
        i1=limlb1**2
        if(magn.eq.2) i1=i1*4
        i2=i1*limlpb
        allocate(lm_tild(nrel_core*limlb1))
        allocate(indxcgr(nrel_core*limlb1,nrel_core*limlb1,2,
     &                   0:3*(magn-1)))
        do ikey=0,1
          if(ikey.eq.0) allocate(jcgr(i2))
          if(ikey.eq.1) allocate(cgr(icg))
          call scgr(icg,ikey)
        enddo
      endif
      allocate(pcor(maxmtcor))
      pcor=0.d0
      allocate(pcor_old(maxmtcor))
      allocate(qcor(maxmtcor))
      allocate(qcor_old(maxmtcor))
      qcor=0.d0
      pcor_old=0.d0
      qcor_old=0.d0
      call me_orb_real
	allocate(n_pbmt0(nsort))
	allocate(lfun_pb(nsort))
	allocate(iopb(natom))
	allocate(iopb0(nsort))
      allocate(gfun(maxwf,nspin))
	allocate(gfund(maxwf,nspin))
      allocate(r(0:maxnrad,nsort))
      allocate(dr(0:maxnrad,nsort))
      allocate(weight(0:maxnrad,nsort))
      allocate(rl(0:maxnrad,0:max_rl,nsort))
      end