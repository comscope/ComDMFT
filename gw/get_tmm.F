      subroutine get_tmm
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: isort,l,iej,j,mtpb,iei,i,mtpbi,ir,ispin,li,je1,le1,li1,
     &	         l1,je3,le3,le,je2,le2,l2,je
	real*8 :: pi8,pi8l,dqdall,c2,deg,fc,ss
      real*8, allocatable ::  work(:),pp(:),tmp(:),f(:),g(:),qq(:),
     &                        pp1(:)
        allocate(work(0:maxnrad),pp(0:maxnrad),tmp(0:maxnrad),
     &	         f(0:maxnrad),g(0:maxnrad))
      c2=clight**2
	pi8=8.d0*pi
	tmm_pb=0.d0
	do isort=1,nsort
c ------ Getting MM-values  ----------------------------------------
        do l=0,lmpb(isort)
	    pi8l=pi8/dfloat(l+l+1)
	    do iej=1,ntle_pb(l,isort)
	        j=ind_prod0(iej,l,isort)
	        mtpb=ind_pb(j,isort)
	        call var_int1(dr(0,isort),fun_pb(mtpb),f,g,tmp(0),
     &	 	                r(0,isort),rl(0,l,isort),work,nrad(isort)+1,
     &                      h(isort))
	      do iei=1,ntle_pb(l,isort)
    	          i=ind_prod0(iei,l,isort)
	          mtpbi=ind_pb(i,isort)
                do ir=0,nrad(isort)
                  pp(ir)=tmp(ir)*fun_pb(mtpbi+ir)*dr(ir,isort)
     &		                                     *r(ir,isort)**2
                enddo
                tmm_pb(iei,iej,l,isort)=pi8l
     &		                           *dqdall(h(isort),pp,nrad(isort))
	      enddo
	    enddo
	    do iej=1,ntle_pb(l,isort)
	      do iei=1,ntle_pb(l,isort)
              tmm_pb(iei,iej,l,isort)=0.5d0*(tmm_pb(iei,iej,l,isort)+
     &                                       tmm_pb(iej,iei,l,isort))
              tmm_pb(iej,iei,l,isort)=tmm_pb(iei,iej,l,isort)
	      enddo
	    enddo
	  enddo
	enddo   !!! over isort
      deallocate(work,pp,tmp,f,g)
c ---------------------------------------------------------
      end