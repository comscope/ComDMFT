      subroutine sigc_im_r(ispin,sig_mi,gr,gc)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: ispin
	real*8, intent(in) :: gr(nr_red,nfun,ndim3_kk(me3_k+1),2)
	complex*16, intent(in) :: gc(nr_red,nfun,ndim3_kk(me3_k+1),2,nrel)
	complex*16, intent(out) :: sig_mi(nr_red,nfun,2,nrel,
     &	                                ndim3_kk(me3_k+1))
	integer :: lm,l,ii,iii,indx,lget,icg1,icg2,icg,jatom,jsort,
     &	         km2,lm2,lf2,km1,lm1,lf1,i,j,ikm2,ir,lf,ie,ind_ir,irl,
     &           ia,isp
	real*8 :: cgg
	complex*16 :: cggc
      isp=min(ispin,nspin_1)
	sig_mi=(0.d0,0.d0)
	do irl=1,nrel
	  do ia=1,2  !! here is over TAU and BETA-TAU
	    do ind_ir=1,ndim3_kk(me3_k+1)
	      ir=n3_mpi_kk(me3_k+1)+ind_ir
	      do jatom=1,natom
	        jsort=is(jatom)
	        do km1=1,lfunm(jsort)
	          lm1=lm_isz(km1,jsort)
	          lf1=lf_isz(km1,jsort)
	          j=io_lem(jatom)+km1-1
	          do km2=1,lfunm(jsort)
	            ikm2=io_lem(jatom)+km2-1
	            lm2=lm_isz(km2,jsort)
	            lf2=lf_isz(km2,jsort)
	            if(irel<=1) then
	              ii = max0(lm1,lm2)
	              iii = min0(lm1,lm2)
                    indx = (ii*(ii-1))/2 + iii
                    icg1 = indxcg(indx)
                    icg2 = indxcg(indx+1) - 1
	            else
		            icg1=indxcgr(lm1,lm2,1,0)
                    icg2=indxcgr(lm1,lm2,2,0)
                  endif  
                  do icg = icg1, icg2
	              if(irel<=1) lm = jcg(icg)
	              if(irel==2) lm = jcgr(icg)
	              l=lget(lm)
	              if(l>lmpb(jsort)) cycle
	              do ie=1,ntle_pb(l,jsort)
	                lf=ind_prod0(ie,l,jsort)
	                i=iopb(jatom)+ind_prod(ie,lm,jsort)-1
	                if(irel<=1) then
	                  cgg=fi0(lf2,lf1,lf,jsort,isp,isp)*cg(icg)
	                  sig_mi(:,j,ia,irl,ind_ir)=
     &	                  sig_mi(:,j,ia,irl,ind_ir)
     &		             -w_rs_mi(:,i,ind_ir)
     &                   *gr(:,ikm2,ind_ir,ia)*cgg
	                else if(irel==2) then
	                  cggc=fi0(lf2,lf1,lf,jsort,isp,isp)*cgr(icg)
	                  sig_mi(:,j,ia,irl,ind_ir)=
     &	                  sig_mi(:,j,ia,irl,ind_ir)
     &			             -w_rs_mi(:,i,ind_ir)
     &                       *gc(:,ikm2,ind_ir,ia,irl)*cggc
                      endif
	              enddo
	            enddo
	          enddo
	        enddo   !! over km1
	      enddo   !! over jatom
	    enddo  !! over ind_k
	  enddo   !! over ia
	enddo   !! over irl
      end