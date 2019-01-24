      subroutine efins(gap)
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer :: n_el,nb,k
	real*8 :: gap,evmax,eemin
      real*8, allocatable :: work(:)
	allocate(work(nbndf*nspin))
      gap=0.d0
	n_el=idnint(nelec)
      if(mod(n_el,2).eq.1.and.nspin.eq.1.and.nrel.eq.1) return
      nb=nelec/2.d0*nrel*nspin
      evmax=-1000.d0
      eemin=1000.d0
      do k=1,npnt
        if(nspin.eq.1) then
          evmax=max(evmax,e_bnd(nb,k,1))
          if(nbndf.gt.nb) eemin=min(eemin,e_bnd(nb+1,k,1))
        else
          call dcopy(nbndf,e_bnd(1,k,1),1,work,1)
          call dcopy(nbndf,e_bnd(1,k,2),1,work(nbndf+1),1)
          call upor(2*nbndf,work)
          evmax=max(evmax,work(nb))
          if(nbndf.gt.nb) eemin=min(eemin,work(nb+1))
        endif
      enddo
      gap=0.d0
      if(nbndf.gt.nb) gap=eemin-evmax
      if(gap.lt.0.d0) gap=0.d0
      if(gap.le.0.d0)return
      if(maswrk.and.iprn.ge.1) then
        write(iun,1010)
        write(iun,1020)gap,gap*13.6058d0
      endif
      efermi=evmax+1.d-06
	deallocate(work)
1010  format(1x,'semiconducting band structure discovered --------->')
1020  format(1x,'dielectric gap value is ',g14.7,' ry ; ',g14.7,' ev')
      end