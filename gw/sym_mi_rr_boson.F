      subroutine sym_mi_rr_boson(vmi,tmp,ind_s)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: ind_s
	real*8, intent(inout) :: vmi(n_pbmt,max_r_star,nqdiv),
     &	                       tmp(n_pbmt,max_r_star,nqdiv)
	integer :: j0,ig,j,j1,iatom,isort,ir,nsta,ir1,jatom,ind0,jnd0,l,
     &	         ie,m,lm,km,lm0
	real*8 :: v(3),v1(3),y0(2*maxpb+1),t0(2*maxpb+1)
	tmp=0.d0
      nsta=num_r_star(ind_s)
	do ir=1,nqdiv
	  do ig=1,ngroup
	    do j=1,nsta
	      j0=ind_r_star(j,ind_s)
	      j1=ip_r(j0,ig)
	      j1=ind_r_star_back(j1)
	      do iatom=1,natom
	        v=rr_coarse(:,ir)+tshift(:,iatom,ig)-tshift_r(:,j0,ig)
	        call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),1)
              call zone1_number(v1,gbas,ndiv,ir1)
	        jatom=ip(iatom,ig)
	        isort=is(iatom)
	        ind0=iopb(iatom)-1   !! zero index for atom iatom
	        jnd0=iopb(jatom)-1
	        do l=0,lmpb(isort)
                do ie=1,ntle_pb(l,isort)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    km=ind_prod(ie,lm,isort)
	              lm0=l+m+1
	              t0(lm0)=vmi(jnd0+km,j1,ir1)
	            enddo
	            call rotate1(t0,y0,l,u(1,ig),1)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    km=ind_prod(ie,lm,isort)
	              lm0=l+m+1
	              tmp(ind0+km,j,ir)=tmp(ind0+km,j,ir)+y0(lm0)
	            enddo
	          enddo  !! over ie
	        enddo   !! over l
	      enddo   !! over iatom
	    enddo  !! over j
	  enddo   !! over ig
	enddo   !! over ir
	vmi=tmp/ngroup
      end