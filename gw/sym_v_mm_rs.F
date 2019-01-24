      subroutine sym_v_mm_rs(ir0,ind_ir0,f_red,f)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ir0,ind_ir0
	complex*16, intent(in) :: f(nqdiv,nlb20,nlb2)
	real*8, intent(out) :: f_red(nlb2,nlb2,maxval(ndim3_k))
	integer :: j,ig,iatom,isort,ig0,ir1,jatom,jsort,jnd0,ind0,ir,
     &	       jnd00,ind00,jnd1,l,m,lm,lm0,i
	real*8 :: y0(2*maxpb+1),t0(2*maxpb+1)
      real*8, allocatable :: tmp(:,:)
      allocate(tmp(limlpb,limlpb))
c -------------------------------------------------------------------
	do jatom=1,natom
	  jsort=is(jatom)
	  jnd0=iob2(jatom)-1
	  do iatom=1,natom
c --------- Find the place of irreducible point ir0 in a common list
          if(ir0>nrr_red(iatom,jatom)) cycle
	    ir=rr_pnt_ind(ir0,iatom,jatom)
	    isort=is(iatom)
	    ind0=iob2(iatom)-1
	    do ig=1,ngroup
	      if(ip(iatom,ig)==iat_1(isort)) then
	        ig0=ig
	        exit
	      endif
	    enddo
		  call r_pair_a_sym(ig0,ir,ir1,iatom,jatom)
	    jnd00=iob2(ip(jatom,ig0))-1
	    ind00=iob20(isort)-1
	    do j=1,(lmpb(jsort)+1)**2
	      jnd1=jnd00+j
	      do l=0,lmpb(isort)
              do m=-l,l
                lm=l*(l+1)+m+1
	          lm0=l+m+1
	          t0(lm0)=f(ir1,ind00+lm,jnd1)
	        enddo
	        call rotate1(t0,y0,l,u(1,ig),1)
              do m=-l,l
                lm=l*(l+1)+m+1
	          lm0=l+m+1
	          tmp(lm,j)=y0(lm0)
	        enddo
	      enddo   !! over l
	    enddo   !! over j
	    do i=1,(lmpb(isort)+1)**2
	      do l=0,lmpb(jsort)
              do m=-l,l
                lm=l*(l+1)+m+1
	          lm0=l+m+1
	          t0(lm0)=tmp(i,lm)
	        enddo
	        call rotate1(t0,y0,l,u(1,ig),1)
              do m=-l,l
                lm=l*(l+1)+m+1
	          lm0=l+m+1
	          f_red(ind0+i,jnd0+lm,ind_ir0)=y0(lm0)
	        enddo
	      enddo   !! over l
	    enddo  !! over i
	  enddo
	enddo
	deallocate(tmp)
      end
