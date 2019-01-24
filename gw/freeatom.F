      subroutine freeatom
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer :: isort,mt,irad
	real*8 :: b,d,ea,rpb
      real*8, allocatable :: rmax(:)
	allocate(nr(nsort))
	allocate(rmax(nsort))
	nrmax=0
	ubi='dft'
      do isort = 1, nsort
c ---- set-up radial mesh with an odd mesh-point at rmt and rmax --
        rmax(isort) =120.d0
	  b=smt(isort)/(exp(h(isort)*nrad(isort))-1.d0)
	  nr(isort)=dlog(rmax(isort)/b+1.d0)/h(isort)
        if (mod(nr(isort),2) .eq. 1) nr(isort) = nr(isort) + 1
        rmax(isort) = b*(dexp(h(isort)*nr(isort)) - 1.d0)
	  nrmax=max(nrmax,nr(isort))
	enddo
      allocate(ro_core(0:nrmax,nsort,nspin_0))
      ro_core=0.d0
      allocate(r_atom(0:nrmax,nsort))
      allocate(ro_atom(0:nrmax,nsort))
c --- We store Info about atomic basis ----------------------------
      do isort=1,nsort
        mt=0
        d=1.d5
	  b=rmax(isort)/(exp(h(isort)*nr(isort))-1.d0)
        ea=dexp(h(isort))
        rpb=b
        do irad=0,nr(isort)
          r_atom(irad,isort)=rpb-b
          rpb=rpb*ea
        enddo
      enddo
c --------------------------------------------------------------------          
	if(iheg) then
	  eny=0.d0
	else
        do isort = 1, nsort
c ---- make atomic sphere selfconsistent ----
	    call atscfp(isort,nr(isort),r_atom(0,isort),ro_atom(0,isort),
     &	              rmax(isort))
	  enddo
	endif
	deallocate(rmax)
      if(nspin==2) then
        do isort=1,nsort
          e_core(:,2,isort)=e_core(:,1,isort)
          eny(:,:,isort,2)=eny(:,:,isort,1)
        enddo
      endif
      end