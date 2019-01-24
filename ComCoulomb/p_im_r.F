      subroutine p_im_r(it0,ind_tau,ispin,p_mi,key,key1,key_p)
c     key = 0 : (G_tau)*(G_beta-tau) contribution to P -----------------
c     key = 1 : (G_tau_UP)*(G_beta-tau_DN) contribution to Xi_T -----
c     key1 = 0 - to use QP Green's function
c     key1 = 1 - to use regular Green's function
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer,intent(in) :: ind_tau,ispin,key,it0,key1,key_p
      double precision, intent(out) ::
     $  p_mi(nr_red,n_pbmt,ndim3_kk(me3_k+1))
      integer :: lm,l,ii,iii,indx,lget,icg1,icg2,icg,jatom,jsort,irl,
     &  km2,lm2,lf2,km1,lm1,lf1,i,ikm1,ikm2,ir,lf,ie,ind_ir,it,
     &  isp,isp1
      double precision :: cgg,dum
      complex*16 :: cggc,dumc
      double precision, allocatable :: gr(:,:,:,:)
      complex*16, allocatable :: gc(:,:,:,:,:)
      if(key==0) then
        isp=min(ispin,nspin_1)
        isp1=isp
      else if(key==1) then
        isp=1
        isp1=nspin
      endif
      if(irel<=1) then
        allocate(gr(nr_red,nfun,ndim3_kk(me3_k+1),2))
        if(key==0) then
          do it=1,2             !! here is over TAU and BETA-TAU
            call g_rs_from_ks_mi(gr,dumc,it,ispin,ind_tau,key1,key_p)
          enddo
        else if(key==1) then
          call g_rs_from_ks_mi(gr,dumc,it0,1,ind_tau,key1,key_p)
          call g_rs_from_ks_mi(gr,dumc,3-it0,nspin,ind_tau,key1,key_p)
        endif
      else if(irel==2) then
        allocate(gc(nr_red,nfun,ndim3_kk(me3_k+1),2,nrel))
        if(key==0) then
          do it=1,2             !! here is over TAU and BETA-TAU
            call g_rs_from_ks_mi(dum,gc,it,ispin,ind_tau,key1,key_p)
          enddo
        endif
      endif
      do ind_ir=1,ndim3_kk(me3_k+1)
        ir=n3_mpi_kk(me3_k+1)+ind_ir
        do jatom=1,natom
          jsort=is(jatom)
          do km1=1,lfunm(jsort)
            lm1=lm_isz(km1,jsort)
            lf1=lf_isz(km1,jsort)
            ikm1=io_lem(jatom)+km1-1
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
                do icg = icg1, icg2
                  lm = jcg(icg)
                  l=lget(lm)
                  if(l>lmpb(jsort)) cycle
                  do ie=1,ntle_pb(l,jsort)
                    lf=ind_prod0(ie,l,jsort)
                    cgg=fi0(lf2,lf1,lf,jsort,isp,isp1)*cg(icg)
                    i=iopb(jatom)+ind_prod(ie,lm,jsort)-1
                    p_mi(:,i,ind_ir)=p_mi(:,i,ind_ir)
     &                -gr(:,ikm1,ind_ir,1)
     &                *gr(:,ikm2,ind_ir,2)*cgg
                  enddo
                enddo
              else if(irel==2) then
                icg1=indxcgr(lm2,lm1,1,0)
                icg2=indxcgr(lm2,lm1,2,0)
                do icg = icg1, icg2
                  lm=jcgr(icg)
                  l=lget(lm)
                  if(l>lmpb(jsort)) cycle
                  do ie=1,ntle_pb(l,jsort)
                    lf=ind_prod0(ie,l,jsort)
                    cggc=fi0(lf2,lf1,lf,jsort,1,1)*cgr(icg)
                    i=iopb(jatom)+ind_prod(ie,lm,jsort)-1
                    do irl=1,2
                      p_mi(:,i,ind_ir)=p_mi(:,i,ind_ir)
     &                  -gc(:,ikm1,ind_ir,1,irl)
     &                  *conjg(gc(:,ikm2,ind_ir,2,irl))
     &                  *cggc
                    enddo
                  enddo
                enddo
              endif
            enddo
          enddo
        enddo
      enddo                     !! over ind_k
      if(irel<=1) deallocate(gr)
      if(irel==2) deallocate(gc)
      end
