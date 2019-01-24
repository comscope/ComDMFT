      subroutine sigc_r_tau_solid
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer :: ispin,ind_k,j,j0,k,i,i0,ia,ib,ic,ind,ind_tau,nbas,
     &               ibas0,irl
      real*8 :: dum
      complex*16 :: dumc
      real*8, allocatable :: gr1(:,:,:,:)
      complex*16, allocatable :: s_c(:,:,:,:),s(:,:),
     &                               s_c1(:,:,:,:,:),gc1(:,:,:,:,:)
c -------- Precalculation of the matrix S*A ----------------------
      allocate(sa(nplw_rm,nbndf,maxval(ndim3_k),nrel,nspin))
      allocate(s(nplw_rm,nbndtm))
      do ind_k=1,ndim3_k(me3_k+1)
        k=n3_mpi_k(me3_k+1)+ind_k
        nbas=nbask(k)/nrel
c ------ We form the overlap matrix -----------------------------
        do j=1,nbas
          j0=indgb(j,k)
          do i=1,nplw_r(k)
            i0=indpw_r(i,k)
            ia=igbs(1,j0)-igbs(1,i0)
            ib=igbs(2,j0)-igbs(2,i0)
            ic=igbs(3,j0)-igbs(3,i0)
            ind=indplw(ia,ib,ic)
            if(complex_ro) then
              s(i,j)=dcmplx(sovr(ind),sovi(ind))
            else
              s(i,j)=dcmplx(sovr(ind),0.d0)
            endif
          enddo
        enddo
        s=s*amega
        do ispin=1,nspin
c -------- Forming S*A ---------------------------------------------
          do irl=1,nrel
            ibas0=(irl-1)*nbas
            call zgemm('n','n',nplw_r(k),n_bnd(k,ispin),nbas,
     &                   (1.d0,0.d0),s,nplw_rm,
     &                 ev_bnd(ibas0+1,1,ind_k,ispin),nbasmpw,
     &                 (0.d0,0.d0),sa(1,1,ind_k,irl,ispin),nplw_rm)
          enddo
        enddo   !! over ispin
      enddo   !! over ind_k
      deallocate(s)
c --------------------------------------------------------------------
      do ind_tau=1,ndim3_tau
          allocate(w_rs_mm(n_pbmtm,n_pbmtm,ndim_nrr(me3_k+1)))
          call w_mm(ind_tau)
          if(ind_tau==1)
     &      call timel('** W_MM_FROM_Q finished ************')
          allocate(w_rs_mi(nr_red,n_pbmt,ndim3_kk(me3_k+1)))
          call v_im_r_from_k(ind_tau)
          if(ind_tau==1)
     &      call timel('** V_IM_R_FROM_K finished **********')
c
c      ---------- MT_MT Part ------------------------------------
c            
          do ispin=1,nspin
            allocate(s_c(nfun,nfun,ndim3_kk(me3_k+1),2))
            s_c=(0.d0,0.d0)
            do ia=1,2   !! here is over TAU and BETA-TAU
              call sigc_mm_r(ispin,ia,ind_tau,s_c(1,1,1,ia))
            enddo
            if(ind_tau==1)
     &        call timel('** SIGC_MM_R finished **************')
            call sigc_mm_k_from_r(ispin,ind_tau,s_c,0)
            if(ind_tau==1)
     &        call timel('** SIGC_MM_K_FROM_R finished *******')
            deallocate(s_c)
c
c      ---------- INT_MT Part ------------------------------------
c      	    
      	    allocate(s_c1(nr_red,nfun,2,nrel,ndim3_kk(me3_k+1)))
            if(irel<=1) then
              allocate(gr1(nr_red,nfun,ndim3_kk(me3_k+1),2))
              do ia=1,2   !! here is over TAU and BETA-TAU
                call g_rs_from_ks_mi(gr1,dumc,ia,ispin,ind_tau)
              enddo
              call sigc_im_r(ispin,s_c1,gr1,dumc)
              deallocate(gr1)
            else if(irel==2) then
              allocate(gc1(nr_red,nfun,ndim3_kk(me3_k+1),2,nrel))
              do ia=1,2   !! here is over TAU and BETA-TAU
                call g_rs_from_ks_mi(dum,gc1,ia,1,ind_tau)
              enddo
              call sigc_im_r(ispin,s_c1,dum,gc1)
              deallocate(gc1)
            endif
            if(ind_tau==1)
     &        call timel('** SIGC_IM_R finished **************')
            call sigc_im_k_from_r(ispin,ind_tau,s_c1)
            if(ind_tau==1)
     &        call timel('** SIGC_IM_K_FROM_R finished *******')
            deallocate(s_c1)
          enddo  !! over ispin
          deallocate(w_rs_mi,w_rs_mm)
c
c      ---------- INT_INT Part ------------------------------------
c      	  
        call sigc_ii_solid(ind_tau)
        if(ind_tau==1)
     &    call timel('** SIGC_II_SOLID finished **********')
      enddo   !! over ind_tau
      deallocate(sa)
      end
