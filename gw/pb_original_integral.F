      subroutine pb_original_integral(ind_k,b)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: ind_k
      complex*16, intent(out) :: b(n_pbtot)
	integer :: isort,iatom,k,iei,i,mtpbi,ir,i3,i2,i1
	real*8 :: pi2,c2,c4,g,sqpi4,q(3),v(3),dqdall,phase,qs
      complex*16 :: cc
	real*8, allocatable :: pp(:)
	k=n3_mpi_k(me3_k+1)+ind_k
      pi2=2.d0*pi
      c2=4.d0*pi/amega
      c4=pi2/par
      sqpi4=sqrt(4.d0*pi)
	b=(0.d0,0.d0)
c ------------------- MT part ------------------------------------
      allocate(pp(0:maxnrad))
      do iatom=1,natom
        isort=is(iatom)	
	  do iei=1,ntle_pb(0,isort)
	    i=ind_prod0(iei,0,isort)
	    mtpbi=ind_pb(i,isort)
          do ir=0,nrad(isort)
            pp(ir)=fun_pb(mtpbi+ir)*dr(ir,isort)*r(ir,isort)**2
          enddo
	    i=iopb(iatom)-1+ind_prod(iei,1,isort)
          b(i)=sqpi4*dqdall(h(isort),pp,nrad(isort))
        enddo
      enddo
      deallocate(pp)
c ------------------- INT part ------------------------------------
      do i=1,nplwgw(k)
        q=pnt(:,k)+gbs(:,indpw_gw(i,k))
        g=c4*sqrt(dot_product(q,q))
        if(g.lt.1.d-6) then
          b(n_pbmt+i)=amegaint
        else
c ------ The whole unit cell contribution -------------------------        
          do i3=0,nrdiv(3)-1
            do i2=0,nrdiv(2)-1
              do i1=0,nrdiv(1)-1
                v=i1*r0b0(:,1)+i2*r0b0(:,2)+i3*r0b0(:,3)
	          phase=pi2*dot_product(q,v)
	          b(n_pbmt+i)=b(n_pbmt+i)+dcmplx(cos(phase),sin(phase))
	        enddo
	      enddo
	    enddo
	    b(n_pbmt+i)=b(n_pbmt+i)*amega/nr_full
c ---------- Subract MT contribution ---------------------------
          do iatom=1,natom
            isort=is(iatom)	
	      phase=pi2*dot_product(q,tau(:,iatom))
	      cc=dcmplx(cos(phase),sin(phase))
	      qs=g*smt(isort)
	      b(n_pbmt+i)=b(n_pbmt+i)-4.d0*pi*cc*(sin(qs)-qs*cos(qs))
     &	                                      /qs**3
          enddo
        endif
      enddo
      end
