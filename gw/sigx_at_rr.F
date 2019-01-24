      subroutine sigx_at_rr(sx,ispin,li0,n)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: ispin,li0,n
      real*8, intent(out) :: sx(0:n,0:n)
	integer :: li,l,i,ntl,ie,lf,je,lf1,ir,l2,jr
	real*8 :: fac,c,pi8,rj
	real*8, allocatable :: cc(:,:)
	pi8=8.d0*pi
	allocate(cc(0:n,0:n))
      sx=0.d0
      do li=1,nrel*lmb(1)+1
	  if(irel.le.1) then
	    l=li-1
	  else
	    l=li/2
	    i=li-2*l
	    if(i.eq.0) i=-1
	  endif
        ntl=ntle(l,1)
        cc=0.d0
        do jr=0,n
          do ie=1,ntl
            lf=indfun0_atom(ie,li)
            c=0.d0
            do je=1,ntl
              lf1=indfun0_atom(je,li)
              c=c+g_beta_atom(ie,je,li,ispin)*pval_atom(jr,lf1,ispin)
            enddo
            do ir=0,n
              cc(ir,jr)=cc(ir,jr)+pval_atom(ir,lf,ispin)*c
            enddo
          enddo
        enddo
	  do l2=0,maxpb
	    fac=lll_at(l2,li0,li)*pi8/(l2+l2+1)
          do jr=1,n
            rj=r_atom(jr,1)**(l2+1)
            do ir=1,jr
              sx(ir,jr)=sx(ir,jr)+fac*cc(ir,jr)*r_atom(ir,1)**l2/rj
            enddo
            rj=r_atom(jr,1)**l2
            do ir=jr+1,n
              sx(ir,jr)=sx(ir,jr)+fac*cc(ir,jr)*rj/r_atom(ir,1)**(l2+1)
            enddo
          enddo
        enddo
      enddo
      deallocate(cc)
      end