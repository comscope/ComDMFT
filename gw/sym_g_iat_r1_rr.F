      subroutine sym_g_iat_r1_rr(f_red,f,iatom,i1,ir)
c     Gives G for given iatom_r'_R 
c     from symmetrized (reduced) representation
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: iatom,i1,ir
	real*8, intent(in) :: f_red(nqdiv,nr_red,nfun)
	real*8, intent(out) :: f(maxel)
	integer :: ig,isort,l,ie,m,lm,km,lm0,in,jn,ind0,i0,ind
	real*8 :: y0(2*maxb+1),t0(2*maxb+1),v(3),v0(3)
      f=0.d0
	isort=is(iatom)
	ind0=io_lem(iatom)-1
	i0=ir_ref(i1)
	ig=r_group(i1)
      v0=rr_pnt(1,ir)*rbas(:,1)+rr_pnt(2,ir)*rbas(:,2)
     &                         +rr_pnt(3,ir)*rbas(:,3)
      v0=v0+tshift(:,iatom,ig)-tshift_r(:,i1,ig)
      call rotate(v0(1),v0(2),v0(3),v(1),v(2),v(3),u(2,ig),1)
      call zone1_number(v,gbas,ndiv,ind)
	do l=0,lmb(isort)
        do ie=1,ntle(l,isort)
	    in=1
	    if(augm(ie,l,isort)/='LOC') in=2
	    do jn=1,in
            do m=-l,l
              lm=l*(l+1)+m+1
              km=indbasa(jn,ie,lm,isort)
	        lm0=l+m+1
	        t0(lm0)=f_red(ind,i0,ind0+km)
	      enddo
	      call rotate1(t0,y0,l,u(1,ig),1)
            do m=-l,l
              lm=l*(l+1)+m+1
              km=indbasa(jn,ie,lm,isort)
	        lm0=l+m+1
	        f(km)=y0(lm0)
	      enddo
	    enddo
	  enddo
	enddo   !! over l
      end