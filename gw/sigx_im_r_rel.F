      subroutine sigx_im_r_rel(ispin,g_mi,sig_mi)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: ispin
	complex*16, intent(in) :: g_mi(nr_red,nfun,ndim3_kk(me3_k+1),nrel)
	complex*16, intent(out) :: sig_mi(nr_red,nfun,nrel,
     &	                                ndim3_kk(me3_k+1))
	integer :: lm,l,lget,icg1,icg2,icg,jatom,jsort,km2,lm2,lf2,km1,
     &	         lm1,lf1,i,j,ikm2,ir,lf,ie,irl,ind_ir,kf1,kf2
	complex*16 :: cgg
	sig_mi=(0.d0,0.d0)
	do irl=1,nrel
	  do ind_ir=1,ndim3_kk(me3_k+1)
	    ir=n3_mpi_kk(me3_k+1)+ind_ir
	    do jatom=1,natom
	      jsort=is(jatom)
	      do km1=1,lfunm(jsort)
	        lm1=lm_isz(km1,jsort)
	        lf1=lf_isz(km1,jsort)
	        kf1=lf1
	        j=io_lem(jatom)+km1-1
	        do km2=1,lfunm(jsort)
	          ikm2=io_lem(jatom)+km2-1
	          lm2=lm_isz(km2,jsort)
	          lf2=lf_isz(km2,jsort)
	          kf2=lf2
		        icg1=indxcgr(lm1,lm2,1,0)
                icg2=indxcgr(lm1,lm2,2,0)
                do icg = icg1, icg2
	            lm = jcgr(icg)
	            l=lget(lm)
	            if(l>lmpb(jsort)) cycle
	            do ie=1,ntle_pb(l,jsort)
	              lf=ind_prod0(ie,l,jsort)
	              cgg=fi0(kf2,kf1,lf,jsort,ispin,ispin)*cgr(icg)
	              i=iopb(jatom)+ind_prod(ie,lm,jsort)-1
	              sig_mi(:,j,irl,ind_ir)=sig_mi(:,j,irl,ind_ir)
     &			                          -s_it_r0(:,i,ind_ir)
     &                                    *g_mi(:,ikm2,ind_ir,irl)
     &                                    *cgg
	            enddo
	          enddo
	        enddo
	      enddo
	    enddo
	  enddo  !! over k
	enddo   !! over irl
      end