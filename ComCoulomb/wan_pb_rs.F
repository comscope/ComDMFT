
      
      subroutine wan2_pb

      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      use wannier90_mod
      implicit none
      include 'mpif.h'
      integer :: iatom,isort,ndim,j,k,ind_k,ispin,kst,jj,
     &  ind,i,n_pbt,i0,isp,jsp,mx2,ind_b,i_pb,
     &  iq,indp,me1,nn,ii,ir,jspin,
     &  jatom,jsort,n1,jnd,ipb,ndim_kk_max,irk,iband1,iband2,iband,
     $  intdum1,intdum2,iwan,ibnd,k0,tempind(1),iorb,jorb,
     $  iopbt,irl,kk,ll,li,mdim,aa,m12,li_i,lim_i,li_j,lim_j,m1,m2,l,
     $  aaindex,ind1,ind2,iw,bb,tempint,
     $  ik,tempcnt1,tempcnt2,tempcnt3,jbnd,kbnd,nj,jpb,iat,jat,n,
     $  ib,jb,korb,iorb_i,jorb_j,ind_i,jnd_j,iq0,jwan
      double precision :: dm,v(3),fc,pi2, tempdble1,tempdble2
      complex*16 :: cc,fac0
      integer, allocatable :: ndim_nr(:),n_mpi_nr(:),tempindmat(:),
     $  wanind_cor(:,:),wanind_iter(:,:),wanind_rest(:,:),
     $  sorting_disen(:,:,:)
      double precision, allocatable :: 
     $  weight_band(:,:,:),tempweightmat(:),pb_int(:,:,:),
     $  tmpmat_dble(:,:,:),eval_disen_sort(:,:,:),
     $  diag(:),tempvec_dble(:)
      complex*16 :: z0(nfun,num_orb_proj,nspin,ndim_kk(me+1)),
     &  a0(nbasmpw,num_orb_proj,nspin,ndim_kk(me+1))


C$$$  $  normalization_check(n_ciw_max,n_ciw_max,2,nspin,natom_ciw)     
      complex*16, allocatable :: t0(:,:),t2(:,:,:),s(:,:,:),tmp(:,:),
     &  zt(:,:,:,:),zi(:,:,:,:,:),b0(:,:),b1(:,:),b2(:,:),
     $  t(:,:),t1(:,:),anew(:,:,:),znew(:,:,:),
     $  fifim(:,:,:,:,:),m_p(:,:,:,:),tempvec(:), hmat_rest(:,:),
     $  tempmat(:,:),cv_mmmm(:,:),cv_mmmm2(:,:),xx_coef(:,:,:),
     $  pi_ww_c(:,:,:,:),tmp2(:),br_mt(:,:,:,:),br_int(:,:,:,:),
     $  tm(:,:,:),tm1(:,:),pi_ww_part(:,:,:,:),
     $  one_m_p(:,:)

      character*50 :: wan_name
!     write(100,*)'n_ciw_atom',n_ciw_atom,'nspin',nspin,'iatom',
!     $ iatom,'n_ciw_atom',n_ciw_atom

!     write(100,*)'n_nu',n_nu,'nporc_tau',nproc_tau,
!     $  'phase_k_rr',phase_k_rr
!     call mpi_bcast(nproc_k,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)


c     write(100,*)'at_orb_ciw',at_orb_ciw,
c     $'enpj_emin',enpj_emin,'enpj_emax',enpj_emax,
c     $'u_mode',u_mode,'p_low_mode',p_low_mode,
c     $'p_local_r',p_local_r,'p_local_k',p_local_k,
c     $'u_linear',u_linear,'projemin',projemin,
c     $'projemax',projemax,'ewald',ewald,'ewald0',ewald0,
c     $'dist_rc',dist_rc,'include_gamma_w',include_gamma_w,
c     $'r_ref',r_ref,'num_orb_cor',num_orb_cor,'num_orb_iter',
c     $num_orb_iter,'num_orb_proj',num_orb_proj,'num_wann',
c     $num_wann,'num_bands',num_bands,'wan_correlated',
c     $wan_correlated,'wan_itinerant',wan_itinerant,
c     $'chem_pot',chem_pot,
c     $'chem_pot_qp',chem_pot_qp,
c     $'chi0_high_k_den',chi0_high_k_den,
c     $'chi0_high_loc',chi0_high_loc,
c     $'chi0_k_den',chi0_k_den,
c     $'chi0_loc',chi0_loc,
c     $'chi0_low_k_den',chi0_low_k_den,
c     $'chi0_low_loc',chi0_low_loc,
c     $'chi_high_k_den',chi_high_k_den,
c     $'chi_high_k_den_tau',chi_high_k_den_tau,
c     $'chi_high_loc',chi_high_loc,
c     $'chi_k_den',chi_k_den,
c     $'chi_k_den_tau',chi_k_den_tau,
c     $'chi_loc',chi_loc,
c     $'ciw_ind',ciw_ind,
c     $'e_qp',e_qp,
c     $'e_qp_dis',e_qp_dis,
c     $'eigenvalues',eigenvalues,
c     $'eigenvalues_cor',eigenvalues_cor,
c     $'eigenvalues_iter',eigenvalues_iter,
c     $'eigenvalues_rest',eigenvalues_rest,
c     $'eps_high_k_den',eps_high_k_den,
c     $'eps_high_loc',eps_high_loc,
c     $'eps_k_den',eps_k_den,
c     $'eps_loc',eps_loc,
c     $'eps_pb_vrt',eps_pb_vrt,
c     $'epsinv_high_k_den',epsinv_high_k_den,
c     $'epsinv_high_loc',epsinv_high_loc,
c     $'epsinv_k_den',epsinv_k_den,
c     $'epsinv_loc',epsinv_loc,
c     $'evolt',evolt,
c     $'gbs',gbs,
c     $'have_disentangled',have_disentangled,
c     $'i_kref',i_kref,
c     $'iatom_ciw',iatom_ciw,
c     $'include_bands',include_bands,
c     $'include_bands_foru',include_bands_foru,
c     $'ind_bands_foru',ind_bands_foru,
c     $'indpw_gw',indpw_gw,
c     $'io_lem',io_lem,
c     $'iopb',iopb,
c     $'irel',irel,
c     $'irvec',irvec,
c     $'is',is,
c     $'iun',iun,
c     $'k_list',k_list,
c     $'k_star',k_star,
c     $'k_star_max',k_star_max,
c     $'lfunm',lfunm,
c     $'lwindow',lwindow,
c     $'m_pi',m_pi,
c     $'magn',magn,
c     $'maswrk',maswrk,
c     $'maxel',maxel,
c     $'me',me,
c     $'me3_k',me3_k,
c     $'me3_tau',me3_tau,
c     $'n3_mpi_k',n3_mpi_k,
c     $'n_bnd',n_bnd,
c     $'n_ciw_atom',n_ciw_atom,
c     $'n_ciw_max',n_ciw_max,
c     $'n_mpi_kk',n_mpi_kk,
c     $'n_pbmt',n_pbmt,
c     $'n_pbmt0',n_pbmt0,
c     $'n_pbmtm',n_pbmtm,
c     $'n_pbtot',n_pbtot,
c     $'natom',natom,
c     $'natom_ciw',natom_ciw,
c     $'nbasmpw',nbasmpw,
c     $'nbndf',nbndf,
c     $'ndegen',ndegen,
c     $'ndim3_k',ndim3_k,
c     $'ndim3_nu',ndim3_nu,
c     $'ndim3_tau',ndim3_tau,
c     $'ndim_kk',ndim_kk,
c     $'ndimwin',ndimwin,
c     $'ndiv',ndiv,
c     $'nfun',nfun,
c     $'nplw_gw',nplw_gw,
c     $'nplwgw',nplwgw,
c     $'npnt',npnt,
c     $'nproc',nproc,
c     $'nqdiv',nqdiv,
c     $'nr_full',nr_full,
c     $'nrdiv',nrdiv,
c     $'nrel',nrel,
c     $'nsort',nsort,
c     $'nspin',nspin,
c     $'nspin_0',nspin_0,
c     $'num_bands',num_bands,
c     $'num_orb_cor',num_orb_cor,
c     $'num_orb_iter',num_orb_iter,
c     $'num_orb_proj',num_orb_proj,
c     $'num_orb_rest',num_orb_rest,
c     $'nwpb',nwpb,
c     $'nwpb_all',nwpb_all,
c     $'nwpb_max',nwpb_max,
c     $'p_low_mode',p_low_mode,
c     $'pi',pi,
c     $'pi_ww',pi_ww,
c     $'pnt',pnt,
c     $'pol_iter_for_u',pol_iter_for_u,
c     $'prod_wan2',prod_wan2,
c     $'q_qp',q_qp,
c     $'q_qp_dis',q_qp_dis,
c     $'r0b0',r0b0,
c     $'r_pnt',r_pnt,
c     $'rb0',rb0,
c     $'rbas',rbas,
c     $'rydberg',rydberg,
c     $'seed_name',seed_name,
c     $'u_k_den',u_k_den,
c     $'u_loc',u_loc,
c     $'u_matrix',u_matrix,
c     $'u_matrix_opt',u_matrix_opt,
c     $'u_mode',u_mode,
c     $'u_static_k',u_static_k,
c     $'v_matrix',v_matrix,
c     $'v_matrix_cor',v_matrix_cor,
c     $'v_matrix_iter',v_matrix_iter,
c     $'v_matrix_rest',v_matrix_rest,
c     $'w_k_den',w_k_den,
c     $'w_loc',w_loc,
c     $'wc_k_den_tau',wc_k_den_tau,
c     $'we_need_u',we_need_u
c     
c     
c     call mpi_bcast(nproc_k,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
c     ------------------------------------------------------------------
!     if(ubi=='dft') then
!     if(allocated(q_qp)) deallocate(q_qp)
!     allocate(q_qp(nbndf,nbndf,ndim3_k(me3_k+1),nspin_0))
!     e_qp=e_bnd
!     z_qp=z_bnd
!     ev_qp=ev_bnd
!     q_qp=(0.d0,0.d0) 
!     chem_pot_qp=chem_pot
!     do ispin=1,nspin
!     do ind_k=1,ndim3_k(me3_k+1)
!     k=n3_mpi_k(me3_k+1)+ind_k
!     do i=1,n_bnd(k,ispin)
!     q_qp(i,i,ind_k,ispin)=dcmplx(1.d0,0.d0) !assuming that calculation is from hartree fock
!     enddo
!     enddo
!     enddo
!     endif

      me1=me+1

      z0=0.0d0
      a0=0.0d0

      if ((magn .eq. 2) .and. (irel .eq. 2)) then
        call all_log
     $    ("wan2_pb is not compatible with irel=2 and magn =2")
        call ending
      endif
      
      if (
     $  (we_need_u) .and.
     $  ((u_mode .eq. 'bnse') .or. (u_mode .eq. 'wnds'))
     $  ) then
        if (allocated(ind_bands_foru)) deallocate(ind_bands_foru)
        allocate(ind_bands_foru(nbndf,npnt,nspin))
        ind_bands_foru=0
        num_orb_rest=num_bands-num_orb_proj
        
      endif
      
      if ((we_need_u) .and. (u_mode .eq. 'wnds')) then
        allocate(eigenvalues_cor(num_orb_cor,nqdiv,nspin))
        eigenvalues_cor=0.0d0
        allocate(v_matrix_cor
     $    (num_bands,num_orb_cor,nqdiv,nspin))
        v_matrix_cor=0.0d0        
        if (num_orb_iter .gt. 0) then
          allocate(eigenvalues_iter(num_orb_iter,nqdiv,nspin))
          eigenvalues_iter=0.0d0          
          allocate(v_matrix_iter
     $      (num_bands,num_orb_iter,nqdiv,nspin))
          v_matrix_iter=0.0d0        
        endif

        allocate(eigenvalues_rest(num_orb_rest,nqdiv,nspin))
        eigenvalues_rest=0.0d0          
        allocate(v_matrix_rest
     $    (num_bands,num_orb_rest,nqdiv,nspin))
        v_matrix_rest=0.0d0

c$$$  allocate(eigenvalues_dis(num_bands,nqdiv,nspin))
c$$$  eigenvalues_dis=0.0d0          
c$$$  allocate(v_matrix_dis
c$$$  $    (num_bands,num_bands,nqdiv,nspin))
c$$$  v_matrix_dis=0.0d0

        allocate(sorting_disen(num_bands,nqdiv,nspin))
        sorting_disen=0
        allocate(eval_disen_sort(num_bands,nqdiv,nspin))
        eval_disen_sort=0.0d0
        allocate(q_qp_dis(nbndf,nbndf,ndim3_k(me3_k+1),nspin_0))
        allocate(e_qp_dis(nbndf,npnt,nspin_0))
        e_qp_dis=e_qp
        q_qp_dis=q_qp
      endif            

      allocate(proj_renorm(num_bands,num_orb_proj,ndim_kk(me+1),nspin))      
      
      do ispin=1, nspin
        if (nspin .eq. 1) then
          seed_name='wannier'
        else
          if (ispin .eq. 1) then
            seed_name='wannier_up'
          else
            seed_name='wannier_dn'
          endif
        endif
        
!     call read_wann_chkpt

#ifdef  USE_HDF5
        call read_wannier_dat_hdf5
#else
        call read_wann_dat
#endif
        
c        write(179,*) 'read_wan_chkpt and bcast done'
c        call flush(179)

!     write(100,*)'eigenvalues',eigenvalues,'lwindow',lwindow,
!     $ 'me',me,'n_mpi_kk',n_mpi_kk,'ndim_kk',ndim_kk,'num_bands',
!     $ num_bands,'num_orb_cor',num_orb_cor,'num_orb_iter',num_orb_iter,
!     $ 'num_orb_proj',num_orb_proj,'num_wann',num_wann,'projemax',
!     $ projemax,'projemin',projemin,'rydberg',rydberg,'u_matrix',
!     $ u_matrix,'u_matrix_opt',u_matrix_opt,'wan_correlated',
!     $ wan_correlated,'wan_itinerant',wan_itinerant

        call projector_renormalization(proj_renorm(1,1,1,ispin))

!     do ind_k=1,ndim_kk(me+1) ! k vector
!     write(100,*) proj_renorm(:,:,ind_k)
!     enddo 

c        write(179,*) 'projector renormalization done'
c        call flush(179)
        do irk = 1, ndim_kk(me+1)        
          call is_orthonormal_dcmplx
     $      (proj_renorm(1,1,irk,ispin),num_bands,
     $      proj_renorm(1,1,irk,ispin),num_bands,
     $      num_orb_proj,num_orb_proj,num_bands,1,
     $      'proj_renorm')
        enddo        

c        write(179,*) 'normalization check'
c        call flush(179)
        call cal_e_imp_rho_bnse(proj_renorm(1,1,1,ispin))                    

        call projector_compare(proj_renorm(1,1,1,ispin))

        allocate(znew(nfun,num_bands,nqdiv)) !
        znew=(0.d0,0.d0)

        allocate(anew(nbasmpw,num_bands,nqdiv)) !
        anew=(0.d0,0.d0)
c        write(179,*)'here-1'
c        call flush(179)

        call za_wan_for_u(ispin,znew,anew)

c        call flush(179)
        

        do iwan=1, num_orb_proj
          do ii=1, nfun
            do ind_k=1,ndim_kk(me+1) ! k vector
              do ibnd=1, num_bands
                k=n_mpi_kk(me+1)+ind_k
                k0=i_kref(k)          
                z0(ii,iwan,ispin,ind_k)
     $            =z0(ii,iwan,ispin,ind_k)
     $            +znew(ii,ibnd,k)*proj_renorm(ibnd,iwan,ind_k,ispin)
              enddo
            enddo
          enddo
        enddo
        deallocate(znew)
        
c        write(179,*)'here1'
c        call flush(179)

        do iwan=1, num_orb_proj
          do ii=1, nbasmpw
            do ind_k=1,ndim_kk(me+1) ! k vector
              do ibnd=1, num_bands
                k=n_mpi_kk(me+1)+ind_k
                k0=i_kref(k)
                a0(ii,iwan,ispin,ind_k)
     $            =a0(ii,iwan,ispin,ind_k)
     $            +anew(ii,ibnd,k)*proj_renorm(ibnd,iwan,ind_k,ispin)
              enddo
            enddo
          enddo
        enddo
        deallocate(anew)

c        write(179,*)'here2'
c        call flush(179)


        if ((we_need_u) .and.
     $    (
     $    (u_mode .eq. 'bnse') .or. (u_mode .eq. 'wnds')        
     $    )
     $    ) then
          

          allocate(weight_band(num_bands,nqdiv,2)) !
          weight_band=0.0d0
          do ind_k=1,ndim_kk(me1)
            ik=n_mpi_kk(me1)+ind_k          
            do ibnd=1, num_bands
              if (num_orb_iter .gt. 0) then
                do iorb=1, num_orb_iter
                  weight_band(ibnd,ik,2)
     $              =weight_band(ibnd,ik,2)
     $              +cdabs(proj_renorm
     $              (ibnd,num_orb_cor+iorb,ind_k,ispin))**2                
                enddo
              endif
              
              do iorb=1, num_orb_cor
                weight_band(ibnd,ik,1)
     $            =weight_band(ibnd,ik,1)
     $            +cdabs(proj_renorm(ibnd,iorb,ind_k,ispin))**2
              enddo
            enddo
          enddo
          if(nproc/=1) call mpi_allreduce_dble
     $      (weight_band,num_bands*nqdiv*2,mpi_sum,MPI_COMM_WORLD)

          if (maswrk) then
            write(iun,*) 'weight projected to correlated_bands'
            do ik=1, npnt
              do ibnd=1, num_bands
                write(iun,'(3i5,2f12.6)') ispin, ik,ibnd,
     $            (e_qp(include_bands(ibnd),ik,ispin)-chem_pot)*
     $            evolt/2.0d0, weight_band(ibnd,ik,1)
              enddo
            enddo
            if (num_orb_iter .gt. 1) then
              write(iun,*) 'weight projected to itinerant_bands'
              do ik=1, npnt
                do ibnd=1, num_bands
                  write(iun,'(3i5,2f12.6)') ispin, ik,ibnd,
     $              (e_qp(include_bands(ibnd),ik,ispin)-chem_pot)*
     $              evolt/2.0d0, weight_band(ibnd,ik,2)
                enddo
              enddo
            endif
          endif

          allocate(tempweightmat(num_bands)) !
          allocate(wanind_cor(num_orb_cor,nqdiv)) !
          wanind_cor=0          
          if (num_orb_iter .gt. 0) then
            allocate(wanind_iter(num_orb_iter,nqdiv)) !
            wanind_iter=0            
          endif
          allocate(wanind_rest(num_orb_rest,nqdiv)) !
          wanind_rest=0
          
          do ik=1, nqdiv
            tempweightmat=weight_band(:,ik,1)
            do ii=1, num_orb_cor
              tempind=maxloc(tempweightmat)
              tempweightmat(tempind(1))=0.0d0
              if ((ik .le. npnt) .and. (u_mode .eq. 'bnse')) then
                ind_bands_foru(include_bands(tempind(1)),ik,ispin)=1
              endif              
              wanind_cor(ii,ik)=tempind(1)
            enddo

            if (num_orb_iter .gt. 0) then

              tempweightmat=weight_band(:,ik,2)

              do ii=1, num_orb_cor
                tempweightmat(wanind_cor(ii,ik))=0.0d0
              enddo
              
              do ii=1, num_orb_iter
                tempind=maxloc(tempweightmat)
                tempweightmat(tempind(1))=0.0d0
                if ((ik .le. npnt) .and. (u_mode .eq. 'bnse')) then                
                  ind_bands_foru(include_bands(tempind(1)),ik,ispin)=2
                endif                
                wanind_iter(ii,ik)=tempind(1)
              enddo
              
            endif

            tempcnt2=0
            do ibnd=1, num_bands
              tempcnt1=0
              do jbnd=1, num_orb_cor
                if (ibnd .eq. wanind_cor(jbnd,ik)) then
                  tempcnt1=tempcnt1+1
                endif
              enddo
              if (num_orb_iter .gt. 0) then
                do jbnd=1, num_orb_iter
                  if (ibnd .eq. wanind_iter(jbnd,ik)) then
                    tempcnt1=tempcnt1+1
                  endif
                enddo
              endif
              if (tempcnt1 .eq. 0) then
                tempcnt2=tempcnt2+1
                wanind_rest(tempcnt2,ik)=ibnd
              endif
            enddo
          enddo

          
          if (num_orb_iter .gt. 0) then
            do ik=1, nqdiv
              do ii=1,num_orb_cor
                do jj=1, num_orb_iter
                  if (wanind_cor(ii,ik)
     $              .eq. wanind_iter(jj,ik)) then
                    if (maswrk) then
                      write(993,*)
     $                  'band_selection for u is wrong, cor, iter',ik
                      write(993,'(1000i5)') wanind_cor(:,ik)
                      write(993,'(1000f12.6)') weight_band(:,ik,1)
                      write(993,'(1000i5)') wanind_iter(:,ik)
                      write(993,'(1000f12.6)') weight_band(:,ik,2)
                    endif
                    call ending
                  endif
                enddo
              enddo
              
              do ii=1,num_orb_rest
                do jj=1, num_orb_iter
                  if (wanind_rest(ii,ik)
     $              .eq. wanind_iter(jj,ik)) then
                    if (maswrk) then
                      write(993,*)
     $                  'band_selection for u is wrong, iter, rest',ik
                      write(993,'(1000i5)') wanind_iter(:,ik)
                      write(993,'(1000f12.6)') weight_band(:,ik,2)
                      write(993,'(1000i5)') wanind_rest(:,ik)                    
                    endif
                    call ending
                  endif
                enddo
              enddo
            enddo
          endif

          do ik=1, nqdiv
            do ii=1,num_orb_rest
              if ((wanind_rest(ii,ik) .le. 0) .or.
     $          (wanind_rest(ii,ik) .gt. num_bands)) then
                if (maswrk) then
                  write(993,*)
     $              'wanind_rest index wrong'
                  write(993,'(1000i5)') wanind_rest(:,ik)                    
                endif
                call ending
              endif
              
              do jj=1, num_orb_cor
                if (wanind_rest(ii,ik)
     $            .eq. wanind_cor(jj,ik)) then
                  if (maswrk) then
                    write(993,*)
     $                'band_selection for u is wrong, cor, rest',ik
                    write(993,'(1000i5)') wanind_cor(:,ik)
                    write(993,'(1000f12.6)') weight_band(:,ik,2)
                    write(993,'(1000i5)') wanind_rest(:,ik)                    
                  endif
                  call ending
                endif
              enddo
            enddo
          enddo
          
          if (u_mode .eq. 'bnse') then
            do ik=1, npnt          
              do ii=1, n_bnd(ik,ispin)
                if (ind_bands_foru(ii,ik,ispin) .eq. 0) then
                  if (e_qp(ii,ik,ispin) .gt.chem_pot_qp) then
                    ind_bands_foru(ii,ik,ispin)=3
                  elseif (e_qp(ii,ik,ispin) .lt.chem_pot_qp) then
                    ind_bands_foru(ii,ik,ispin)=4
                  endif
                endif
              enddo            
            enddo

          endif          
          

          if (
     $      (we_need_u) .and.
     $      ((u_mode .eq. 'bnse') .or. (u_mode .eq. 'wnds'))          
     $      ) then          
            
            if (maswrk) then
              write(iun,*) 'correlated bands from band selection'
              do ik=1, nqdiv
                write(iun,'(i5, 4x,100(i5,f12.6,2x))')
     $            ik, (wanind_cor(ii, ik),
     $            eigenvalues(wanind_cor(ii, ik), ik),
     $            ii=1, num_orb_cor)
              enddo
              write(iun,*) 'itinerant bands from band selection'
              do ik=1, nqdiv
                write(iun,'(i5, 4x,100(i5,f12.6,2x))')
     $            ik, (wanind_iter(ii, ik),
     $            eigenvalues(wanind_iter(ii, ik), ik),
     $            ii=1, num_orb_iter)
              enddo              
              write(iun,*) 'the rest bands from band selection'
              do ik=1, nqdiv
                write(iun,'(i5, 4x,100(i5,f12.6,2x))')
     $            ik, (wanind_rest(ii, ik),
     $            eigenvalues(wanind_rest(ii, ik), ik),
     $            ii=1, num_orb_rest)
              enddo                            
            endif
          endif
          

          if ((we_need_u) .and. (u_mode .eq. 'wnds')) then
            
            allocate(tempvec(num_bands)) !


            do ind_k=1, ndim_kk(me+1)
              ik=n_mpi_kk(me+1)+ind_k
              allocate(tempmat(num_bands,num_bands)) !
              tempmat=0.0d0
              call wan_disentangled_hamiltonian
     $          (num_orb_cor,num_bands,
     $          eigenvalues(1,ik),proj_renorm(1,1,ind_k,ispin),
     $          eigenvalues_cor(1,ik,ispin),v_matrix_cor(1,1,ik,ispin))
              do ibnd=1, num_orb_cor
                do jbnd=1, num_bands
                  tempmat(jbnd,ibnd)
     $              =v_matrix_cor(jbnd,ibnd,ik,ispin)
                enddo
              enddo
              

              if (num_orb_iter .gt. 0) then            
                call wan_disentangled_hamiltonian
     $            (num_orb_iter,num_bands,
     $            eigenvalues(1,ik),
     $            proj_renorm(1,(num_orb_cor+1),ind_k,ispin),
     $            eigenvalues_iter(1,ik,ispin),
     $            v_matrix_iter(1,1,ik,ispin))
                do ibnd=1, num_orb_iter
                  do jbnd=1, num_bands
                    tempmat(jbnd,ibnd+num_orb_cor)
     $                =v_matrix_iter(jbnd,ibnd,ik,ispin)
                  enddo
                enddo                
              endif

              allocate(one_m_p(num_bands,num_bands)) !
              one_m_p=0.0d0

              allocate(tempvec_dble(num_bands))
              tempvec_dble=0.0d0
              do ibnd=1, num_bands
                one_m_p(ibnd,ibnd)=1.0d0
              enddo
              
! projection operator to the wannier space
              call zgemm('n','c',num_bands,num_bands,
     $          num_orb_proj,(-1.d0,0.d0),
     $          tempmat,num_bands,
     $          tempmat,num_bands,(1.d0,0.d0),
     $          one_m_p,num_bands)
              
              call hermitianeigen_cmplxdouble
     $          (num_bands,tempvec_dble,one_m_p)
c              write(179,'(a, i5, 10000f12.6)')
c     $          'one_m_p eigenvalues', ik, tempvec_dble(1:num_orb_proj)
c              write(179,'(a, i5, 10000f12.6)')
c     $          'one_m_p eigenvalues', ik,
c     $          tempvec_dble((num_orb_proj+1):num_bands)
              do ibnd=1, num_orb_rest
                tempmat(:,ibnd+num_orb_proj)
     $            =one_m_p(:,ibnd+num_orb_proj)
              enddo
c$$$              
c$$$              do ibnd=1, num_orb_rest
c$$$                tempmat(wanind_rest(ibnd,ik),ibnd+num_orb_proj)
c$$$     $            =1.0d0
c$$$              enddo  
c$$$!     modified gram-schmidt 
c$$$              
c$$$              do ibnd=1, num_orb_rest
c$$$                tempvec
c$$$     $            =tempmat(:,ibnd+num_orb_proj)
c$$$                do jbnd=1,num_orb_proj+ibnd-1
c$$$                  tempvec=
c$$$     $              tempvec
c$$$     $              -sum(dconjg(tempmat(:,jbnd))
c$$$     $              *tempvec)
c$$$     $              *tempmat(:,jbnd)
c$$$                enddo
c$$$                tempmat(:,ibnd+num_orb_proj)
c$$$     $            =tempvec/dsqrt(sum(cdabs(tempvec)**2))
c$$$              enddo

              deallocate(one_m_p)
              deallocate(tempvec_dble)
!     orthonormalization check
              call is_orthonormal_dcmplx
     $          (tempmat,num_bands,
     $          tempmat,num_bands,
     $          num_bands,num_bands,num_bands,1,
     $          'basis_after_orthonormal_basis_set_construction')

              call all_log('gram-schmidt done')

!     disentangle the rest bands
              allocate(hmat_rest(num_orb_rest, num_orb_rest)) !
              hmat_rest=0.0d0
              do ibnd=1, num_orb_rest
                do jbnd=1, num_orb_rest
                  do kbnd=1, num_bands
                    hmat_rest(ibnd,jbnd)=hmat_rest(ibnd,jbnd)
     $                +dconjg(tempmat
     $                (kbnd,ibnd+num_orb_proj))
     $                *eigenvalues(kbnd,ik)
     $                *tempmat
     $                (kbnd,jbnd+num_orb_proj)
                  enddo
                enddo
              enddo

              call hermitianeigen_cmplxdouble(num_orb_rest,
     $          eigenvalues_rest(:,ik,ispin),hmat_rest)

              call zgemm('n','n',num_bands,num_orb_rest,num_orb_rest,
     $          (1.d0,0.d0),tempmat(1,1+num_orb_proj),num_bands,
     &          hmat_rest,num_orb_rest,(0.d0,0.d0),
     $          v_matrix_rest(1,1,ik,ispin),num_bands)

              call all_log('h_rest done')

              deallocate(hmat_rest)
              
              deallocate(tempmat)

              do iwan=1, num_orb_cor
                eval_disen_sort(iwan,ik,ispin)
     $            =eigenvalues_cor(iwan,ik,ispin)
              enddo
              do iwan=1, num_orb_iter
                eval_disen_sort(iwan+num_orb_cor,ik,ispin)
     $            =eigenvalues_iter(iwan,ik,ispin)
              enddo
              do iwan=1, num_orb_rest
                eval_disen_sort
     $            (iwan+num_orb_proj,ik,ispin)
     $            =eigenvalues_rest(iwan,ik,ispin)
              enddo
              
              call all_log('energy_sorting start')
              call sort_dble(eval_disen_sort(1,ik,ispin),
     $          num_bands,sorting_disen(1,ik,ispin))
              
              call all_log('energy_sorting done')

!     orthonormalization check

            enddo
            call all_log('disentangled hamiltonian constructed')

            if(nproc/=1) then
              call mpi_allreduce_dble
     $          (eigenvalues_cor(:,:,ispin),num_orb_cor*nqdiv,
     $          mpi_sum,MPI_COMM_WORLD)
              if(maswrk) then
                write(iun,*) 'eigenvalues_cor'
                do ik=1, nqdiv
                  write(iun,'(i5, 100f12.6)')
     $              ik, (eigenvalues_cor(iwan,ik,ispin),
     $              iwan=1, num_orb_cor)
                enddo
              endif
              call mpi_allreduce_dcmplx
     $          (v_matrix_cor(:,:,:,ispin),num_orb_cor*num_bands*nqdiv,
     $          mpi_sum,MPI_COMM_WORLD)
              if (num_orb_iter .gt. 0) then
                call mpi_allreduce_dble
     $            (eigenvalues_iter(:,:,ispin),num_orb_iter*nqdiv,
     $            mpi_sum,MPI_COMM_WORLD)
                if(maswrk) then
                  write(iun,*) 'eigenvalues_iter'
                  do ik=1, nqdiv
                    write(iun,'(i5, 100f12.6)')
     $                ik, (eigenvalues_iter(iwan,ik,ispin),
     $                iwan=1, num_orb_iter)
                  enddo                  
                endif                
                call mpi_allreduce_dcmplx
     $            (v_matrix_iter(:,:,:,ispin),
     $            num_orb_iter*num_bands*nqdiv,
     $            mpi_sum,MPI_COMM_WORLD)
              endif
              call mpi_allreduce_dble
     $          (eigenvalues_rest(:,:,ispin),num_orb_rest*nqdiv,
     $          mpi_sum,MPI_COMM_WORLD)
              if(maswrk) then
                write(iun,*) 'eigenvalues_rest'
                do ik=1, nqdiv
                  write(iun,'(i5, 100f12.6)')
     $              ik, (eigenvalues_rest(iwan,ik,ispin),
     $              iwan=1, num_orb_rest)
                enddo                                  
              endif              
              call mpi_allreduce_dcmplx
     $          (v_matrix_rest(:,:,:,ispin),
     $          num_orb_rest*num_bands*nqdiv,
     $          mpi_sum,MPI_COMM_WORLD)

              call mpi_allreduce_int
     $          (sorting_disen(:,:,ispin),
     $          num_bands*nqdiv,
     $          mpi_sum,MPI_COMM_WORLD)
              call mpi_allreduce_dble
     $          (eval_disen_sort(:,:,ispin),
     $          num_bands*nqdiv,
     $          mpi_sum,MPI_COMM_WORLD)                            
            endif
!     construction of v_matrix_dis and eigenvalues_dis




c$$$  call wannier_interpolation_path_disentangle(proj_renorm,1)
c$$$  if (num_orb_iter .gt. 0) then
c$$$  call wannier_interpolation_path_disentangle(proj_renorm,2)
c$$$  endif
c$$$  
c$$$  call wannier_interpolation_grid_pdos_disentangle
c$$$  $        (proj_renorm,1)
c$$$  
c$$$  if (num_orb_iter .gt. 0) then            
c$$$  call wannier_interpolation_grid_pdos_disentangle
c$$$  $          (proj_renorm,2)
c$$$  endif

            

            if (maswrk) then
              do ik=1, nqdiv
                call is_orthonormal_dcmplx
     $            (v_matrix_cor(1,1,ik,ispin),num_bands,
     $            v_matrix_cor(1,1,ik,ispin),num_bands,
     $            num_orb_cor,num_orb_cor,num_bands,1,'cor_cor')
                if (num_orb_iter.gt.0) then
                  call is_orthonormal_dcmplx
     $              (v_matrix_cor(1,1,ik,ispin),num_bands,
     $              v_matrix_iter(1,1,ik,ispin),num_bands,
     $              num_orb_cor,num_orb_iter,num_bands,0,'cor_iter')
                  call is_orthonormal_dcmplx
     $              (v_matrix_iter(1,1,ik,ispin),num_bands,
     $              v_matrix_iter(1,1,ik,ispin),num_bands,
     $              num_orb_iter,num_orb_iter,num_bands,1,'iter_iter')
                  call is_orthonormal_dcmplx
     $              (v_matrix_iter(1,1,ik,ispin),num_bands,
     $              v_matrix_rest(1,1,ik,ispin),num_bands,
     $              num_orb_iter,num_orb_rest,num_bands,0,'iter_rest')
                endif
                call is_orthonormal_dcmplx
     $            (v_matrix_cor(1,1,ik,ispin),num_bands,
     $            v_matrix_rest(1,1,ik,ispin),num_bands,
     $            num_orb_cor,num_orb_rest,num_bands,0,'cor_rest')

                call is_orthonormal_dcmplx
     $            (v_matrix_rest(1,1,ik,ispin),num_bands,
     $            v_matrix_rest(1,1,ik,ispin),num_bands,
     $            num_orb_rest,num_orb_rest,num_bands,1,'rest_rest')                                           
              enddo
            endif

            call all_log('orthonormaliztion check of v_matrix')
!     write(100,*)'npnt',npnt,'num_orb_cor',num_orb_cor,
!     $ 'num_orb_iter',num_orb_iter,'num_orb_rest',num_orb_rest,
!     $ 'chem_pot',chem_pot,'num_bands',num_bands,'nqdiv',nqdiv,
!     $ 'nspin',nspin,'sorting_disen',sorting_disen,'eigenvalues',
!     $ eigenvalues,'proj_renorm',proj_renorm,'eigenvalues_cor',
!     $ eigenvalues_cor,'v_matrix_cor',v_matrix_cor,'eigenvalues_iter',
!     $ eigenvalues_iter,'v_matrix_iter',v_matrix_iter,
!     $ 'eigenvalues_rest',eigenvalues_rest
            do ik=1, npnt
              do iwan=1, num_orb_cor
                e_qp_dis(include_bands(
     $            sorting_disen(iwan,ik,ispin)),ik,ispin)
     $            =eigenvalues_cor(iwan,ik,ispin)/rydberg+chem_pot
              enddo
              if (num_orb_iter .gt. 0) then
                do iwan=1, num_orb_iter
                  e_qp_dis(include_bands(
     $              sorting_disen(iwan+num_orb_cor,ik,ispin)),
     $              ik,ispin)
     $              =eigenvalues_iter(iwan,ik,ispin)/rydberg+chem_pot
                enddo
              endif
              do iwan=1, num_orb_rest
                e_qp_dis(include_bands(
     $            sorting_disen(iwan+num_orb_proj,ik,ispin)),ik,ispin)
     $            =eigenvalues_rest(iwan,ik,ispin)/rydberg+chem_pot
              enddo
            enddo


            call all_log('e_qp_dis construction')
            
            allocate(tempmat(nbndf,num_bands))
            
            do ind_k=1,ndim3_k(me3_k+1) 
              ik=n3_mpi_k(me3_k+1)+ind_k 
              tempmat=0.0d0
              call zgemm('n','n',nbndf,num_orb_cor,num_bands,
     $          (1.d0,0.d0),q_qp(1,include_bands(1),ind_k,ispin),nbndf,
     &          v_matrix_cor(1,1,ik,ispin),num_bands,(0.d0,0.d0),
     $          tempmat,nbndf)
              if (num_orb_iter .gt. 0) then
                call zgemm('n','n',nbndf,num_orb_iter,num_bands,
     $            (1.d0,0.d0),q_qp(1,include_bands(1),ind_k,ispin),
     $            nbndf,v_matrix_iter(1,1,ik,ispin),num_bands,
     $            (1.d0,0.d0),tempmat(1,num_orb_cor+1),nbndf)
              endif
              call zgemm('n','n',nbndf,num_orb_rest,num_bands,
     $          (1.d0,0.d0),q_qp(1,include_bands(1),ind_k,ispin),nbndf,
     &          v_matrix_rest(1,1,ik,ispin),num_bands,(1.d0,0.d0),
     $          tempmat(1,num_orb_proj+1),nbndf)
              do ibnd=1, num_bands
                do jbnd=1, n_bnd(ik,ispin)
                  q_qp_dis(jbnd,
     $              include_bands(
     $              sorting_disen(ibnd,ik,ispin)),
     $              ind_k,ispin)=tempmat(jbnd,ibnd)
                enddo
              enddo
            enddo
            deallocate(tempmat)

            call all_log('q_qp_dis construction')            

            if (maswrk) then
              write(iun,*) 'e_qp comparizon after disentanblement'
              do ik=1, npnt
                do ibnd=1, n_bnd(ik,ispin)
                  write(iun,'(3i5, 3f12.6)') ispin,ik,ibnd,
     $              (e_qp_dis(ibnd,ik,ispin)-chem_pot)*rydberg,
     $              (e_qp(ibnd,ik,ispin)-chem_pot)*rydberg,
     $              (e_qp_dis(ibnd,ik,ispin)-chem_pot)*rydberg-
     $              (e_qp(ibnd,ik,ispin)-chem_pot)*rydberg                  
                enddo
              enddo
            endif

            do ind_k=1,ndim3_k(me3_k+1) 
              ik=n3_mpi_k(me3_k+1)+ind_k             
              tempint=include_bands(1)-1
              call is_orthonormal_dcmplx
     $          (q_qp_dis(1,1,ind_k,ispin),nbndf,
     $          q_qp(1,1,ind_k,ispin),nbndf,
     $          tempint,tempint,n_bnd(ik,ispin),1,
     $          'q_qp_dis_low, q_qp_low')
              tempcnt1=n_bnd(ik,ispin)-include_bands(num_bands)
              tempcnt2=include_bands(num_bands)+1              
              call is_orthonormal_dcmplx
     $          (q_qp_dis(1,tempcnt2,ind_k,ispin),nbndf,
     $          q_qp(1,tempcnt2,ind_k,ispin),nbndf,
     $          tempcnt1,tempcnt1,n_bnd(ik,ispin),1,
     $          'q_qp_dis_high, q_qp_high')
              tempcnt1=n_bnd(ik,ispin)-include_bands(num_bands)
              tempcnt2=include_bands(num_bands)+1
              tempcnt3=include_bands(1)-1
              call is_orthonormal_dcmplx
     $          (q_qp_dis(1,tempcnt2,ind_k,ispin),nbndf,
     $          q_qp(1,1,ind_k,ispin),nbndf,
     $          tempcnt1,tempcnt3,n_bnd(ik,ispin),0,
     $          'q_qp_dis_high, q_qp_low')
              call is_orthonormal_dcmplx
     $          (q_qp(1,tempcnt2,ind_k,ispin),nbndf,
     $          q_qp_dis(1,1,ind_k,ispin),nbndf,
     $          tempcnt1,tempcnt3,n_bnd(ik,ispin),0,
     $          'q_qp_dis_low, q_qp_high')
              tempcnt2=include_bands(num_bands)+1
              call is_orthonormal_dcmplx
     $          (q_qp_dis(1,tempcnt2,ind_k,ispin),nbndf,
     $          q_qp_dis(1,tempcnt2,ind_k,ispin),nbndf,
     $          num_bands,num_bands,n_bnd(ik,ispin),1,
     $          'q_qp_dis_wan, q_qp_dis_wan')                            
            enddo

            call all_log('q_qp_dis orthonormalization check done')

            do ik=1, npnt
              do ii=1, num_orb_cor
                ind_bands_foru(include_bands(
     $            sorting_disen(ii,ik,ispin)),
     $            ik,ispin)=1
              enddo
            enddo
            
            if (num_orb_iter .gt. 0) then
              do ik=1, npnt
                do ii=1, num_orb_iter
                  ind_bands_foru(
     $              include_bands(
     $              sorting_disen(ii+num_orb_cor,ik,ispin)),
     $              ik,ispin)=2
                enddo
              enddo                
            endif
            do ik=1, npnt          
              do ii=1, n_bnd(ik,ispin)
                if (ind_bands_foru(ii,ik,ispin) .eq. 0) then
                  if (e_qp_dis(ii,ik,ispin) .gt.chem_pot_qp) then
                    ind_bands_foru(ii,ik,ispin)=3
                  elseif (e_qp_dis(ii,ik,ispin) .lt.chem_pot_qp) then
                    ind_bands_foru(ii,ik,ispin)=4
                  endif
                endif
              enddo            
            enddo
            call all_log('ind_bands_foru construction')
            deallocate(tempvec)
          endif  
          
          deallocate(tempweightmat)
          deallocate(wanind_cor)
          if (num_orb_iter .gt. 0) then          
            deallocate(wanind_iter)
          endif
          deallocate(wanind_rest)
          deallocate(weight_band)

        endif
      enddo
      
      if (maswrk) then
        write(iun,'(a)') 'ind_bands_foru'
        if (u_mode .eq. 'bnse') then
          do ispin=1, nspin
            do ik=1,npnt
              do ii=1, n_bnd(ik,ispin)
                write(iun,'(4i5,2f12.6)') ispin,ik,ii,
     $            ind_bands_foru(ii,ik,ispin),
     $            (e_qp(ii,ik,ispin)-chem_pot)*evolt/2.0d0
              enddo
            enddo
          enddo
          write(iun,'(a)') 'ind_band_corr,ispin,ik,index, energy'
          do ispin=1,nspin          
            do ik=1,npnt
              do ii=1,n_bnd(ik,ispin)
                if (ind_bands_foru(ii,ik,ispin) .eq. 1) then
                  write(iun,'(3i5, 2f12.6)') ispin,ik,ii,
     $              (e_qp(ii,ik,ispin)-chem_pot)*evolt/2.0d0
                endif
              enddo
            enddo
          enddo
          
          if (num_orb_iter .gt. 0) then                            
            write(iun,'(a)') 'ind_band_iter,ispin,ik,index, energy'
            do ispin=1,nspin
              do ik=1,npnt
                do ii=1,n_bnd(ik,ispin)
                  if (ind_bands_foru(ii,ik,ispin) .eq. 2) then
                    write(iun,'(3i5, 2f12.6)') ispin,ik,ii,
     $                (e_qp(ii,ik,ispin)-chem_pot)*evolt/2.0d0
                  endif
                enddo
              enddo
            enddo
          endif                  
        elseif (u_mode .eq. 'wnds') then
          do ispin=1, nspin
            do ik=1,npnt
              do ii=1, n_bnd(ik,ispin)
                write(iun,'(4i5,2f12.6)') ispin,ik,ii,
     $            ind_bands_foru(ii,ik,ispin),
     $            (e_qp(ii,ik,ispin)-chem_pot)*evolt/2.0d0,
     $            (e_qp_dis(ii,ik,ispin)-chem_pot)*evolt/2.0d0              
              enddo
            enddo
          enddo
          write(iun,'(a)') 'ind_band_corr,ispin,ik,index, energy'
          do ispin=1,nspin          
            do ik=1,npnt
              do ii=1,n_bnd(ik,ispin)
                if (ind_bands_foru(ii,ik,ispin) .eq. 1) then
                  write(iun,'(3i5, 2f12.6)') ispin,ik,ii,
     $              (e_qp(ii,ik,ispin)-chem_pot)*evolt/2.0d0,
     $              (e_qp_dis(ii,ik,ispin)-chem_pot)*evolt/2.0d0                
                endif
              enddo
            enddo
          enddo
          
          if (num_orb_iter .gt. 0) then                            
            write(iun,'(a)') 'ind_band_iter,ispin,ik,index, energy'
            do ispin=1,nspin
              do ik=1,npnt
                do ii=1,n_bnd(ik,ispin)
                  if (ind_bands_foru(ii,ik,ispin) .eq. 2) then
                    write(iun,'(3i5, 2f12.6)') ispin,ik,ii,
     $                (e_qp(ii,ik,ispin)-chem_pot)*evolt/2.0d0,
     $                (e_qp_dis(ii,ik,ispin)-chem_pot)*evolt/2.0d0                  
                  endif
                enddo
              enddo
            enddo
          endif          
        endif
        
      endif

      
!!!!!!!!!!!!!!!!start here
      fac0=(1.d0,0.d0)/dble(nqdiv)
      pi2=pi+pi
      allocate(ndim_nr(nproc))
      allocate(n_mpi_nr(nproc))
      call size_shift_par(nr_full,nproc,ndim_nr,n_mpi_nr)
      mx2=n_ciw_max**2
      allocate(fifim(maxel,maxel,n_pbmtm,nsort,nspin))
      allocate(pb_int(n_pbmtm,nsort,nspin))
c     ---------- Interaction matrices -----------------------------------
      call fi0_full_new(fifim)
c     -------- Overlap between bi-products of local Wannier functions -

      allocate(zt(nfun,num_orb_proj,nspin,ndim_kk(me1)))      
      call z_wan_rr_mt(z0,zt)
      allocate(zi(nr_full,num_orb_proj,nrel,nspin,ndim_kk(me1)))      
      call z_wan_rr_int(a0,zi)
      
!     write(100,*)'a0',a0,'zi',zi,'amega',amega,'gbs',gbs,
!     $ 'i_kref',i_kref,'indgb',indgb,'n_mpi_kk',n_mpi_kk,
!     $ 'nbask',nbask,'nbasmpw',nbasmpw,'ndim_kk',ndim_kk,
!     $ 'ndiv',ndiv,'nproc',nproc,'nqdiv',nqdiv,'nr_full',
!     $ nr_full,'nrdiv',nrdiv,'nrel',nrel,'nspin',nspin,
!     $ 'num_orb_proj',num_orb_proj,'pnt',pnt,'r0b0',r0b0,
!     $ 'r_pnt',r_pnt,'rb0',rb0,'rbas',rbas,'nfun',nfun,
!     $ 'cg',cg,'cgr',cgr,'fi0',fi0,'ind_prod',ind_prod,
!     $ 'ind_prod0',ind_prod0,'indxcg',indxcg,'indxcgr',
!     $ indxcgr,'irel',irel,'jcg',jcg,'jcgr',jcgr,'lf_isz',
!     $ lf_isz,'lfunm',lfunm,'lm_isz',lm_isz,'lmpb',lmpb,
!     $ 'maxel',maxel,'n_pbmtm',n_pbmtm,'nsort',nsort,
!     $ 'nspin',nspin,'ntle_pb',ntle_pb
!     call mpi_bcast(irel,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
c     ------ We form the overlap matrix -----------------------------

      allocate(br_mt(n_pbmt,n_ciw_max**2*nspin,ndim_kk(me1),
     $  natom_ciw))
      br_mt=0.0d0

      allocate(br_int(nplw_gw,n_ciw_max**2*nspin,ndim_kk(me1),
     $  natom_ciw))
      br_int=0.0d0

      allocate(s(nplw_gw,nplw_gw,ndim_kk(me1)))


      allocate(b0(nqdiv,nr_full))
      allocate(b1(nr_full,nqdiv))
      allocate(b2(nr_full,nqdiv))      


c     --------------- INT contribution -------------------------------

      allocate(prod_wan2(n_pbtot,num_orb_proj,num_orb_proj,
     $  nspin,ndim3_k(me3_k+1),k_star_max))
      prod_wan2=0.0d0

!     write(100,*)'at_orb_ciw',at_orb_ciw,'num_orb_proj',
!     $ num_orb_proj,'nr_full',nr_full,'zi',zi,'nqdiv',nqdiv,
!     $ 'ndim_nr',ndim_nr,'n_mpi_nr',n_mpi_nr,'r_pnt',r_pnt,
!     $ 'r0b0',r0b0,'ndiv',ndiv
!     call mpi_bcast(nproc_k,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)     
      
      do ispin=1,nspin
        do iorb=1,num_orb_proj  !2
          iat=at_orb_ciw(1, iorb)
          iorb_i=at_orb_ciw(2, iorb)
          do jorb=1,num_orb_proj !3
            jat=at_orb_ciw(1, jorb)
            jorb_j=at_orb_ciw(2, jorb)

            b0=(0.d0,0.d0)
            do ii=1,nr_full
              do irl=1, nrel
                do ind_k=1,ndim_kk(me1)
                  ir=n_mpi_kk(me1)+ind_k
                  b0(ir,ii)
     $              =dconjg(zi(ii,iorb,irl,ispin,ind_k))
     $              *zi(ii,jorb,irl,ispin,ind_k)
                enddo
              enddo
            enddo
            if(nproc/=1) call mpi_allreduce_dcmplx
     $        (b0,nqdiv*nr_full,mpi_sum,MPI_COMM_WORLD)
            b1=(0.d0,0.d0)
            do ind_i=1,ndim_nr(me1)
              ii=n_mpi_nr(me1)+ind_i
              v(:)=r_pnt(1,ii)*r0b0(:,1)+r_pnt(2,ii)*r0b0(:,2)
     &          +r_pnt(3,ii)*r0b0(:,3)
              call fft3(ndiv(1),ndiv(2),ndiv(3),1,b0(1,ii),-1)
              do iq=1,nqdiv
                call zone1_number(pnt(:,iq),rb0,ndiv,ind)
c     ----- Multiplication with the phase factors ------------------
                fc=pi2*dot_product(pnt(:,iq),v)
                cc=cdexp(dcmplx(0.d0,-fc))
                b1(ii,iq)=b0(ind,ii)*cc
              enddo
            enddo
            if(nproc/=1) call mpi_allreduce_dcmplx
     $        (b1,nqdiv*nr_full,mpi_sum,MPI_COMM_WORLD)
            
            b2=(0.d0,0.d0)
            do ind_k=1,ndim_kk(me1)
              iq=n_mpi_kk(me1)+ind_k
              k0=i_kref(iq)                        
              
              call fft3(nrdiv(1),nrdiv(2),nrdiv(3),1,b1(1,iq),-1)
              b2(:,iq)=b1(:,iq)/dble(nr_full)
              
              if (iat .eq. jat) then
                iatom=iatom_ciw(iat)
                isort=is(iatom)
                n=n_ciw_atom(iatom)
                nn=n**2                
                call indexing(nn*nspin,3,(/n,n,nspin/),1,korb,
     $            (/iorb_i,jorb_j,ispin/))
                do ii=1,nplwgw(k0)
                  i0=indpw_gw(ii,iq)
                  call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
                  br_int(ii,korb,ind_k,iat)=b2(ind,iq)
                enddo
              endif
            enddo
            if(nproc/=1) call mpi_allreduce_dcmplx
     $        (b2,nqdiv*nr_full,mpi_sum,MPI_COMM_WORLD)

            do ind_k=1,ndim3_k(me3_k+1) 
              k=n3_mpi_k(me3_k+1)+ind_k 
              do kst=1,k_star(k)
                iq=k_list(kst,k)
                do ii=1,nplwgw(k)
                  i0=indpw_gw(ii,iq)
                  call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
                  prod_wan2(n_pbmt+ii,iorb,jorb,ispin,ind_k,kst)
     $              =b2(ind,iq)
                enddo
              enddo
            enddo
          enddo                 !! over i
        enddo                   !! over j
      enddo                     !! over ispin

c     --------------- MT contribution -------------------------------
      allocate(t(num_orb_proj,maxel))
      allocate(t1(num_orb_proj,num_orb_proj))
      allocate(t2(nqdiv,num_orb_proj,num_orb_proj))
      t2=0.0d0
      
      do ispin=1,nspin      
        do jatom=1,natom
          jsort=is(jatom)
          n1=lfunm(jsort)
          jnd=io_lem(jatom)     !! first index in NFUN
          do jj=1,n_pbmt0(jsort)
            jpb=iopb(jatom)+jj-1
            
            t2=0.0d0
            do ind_k=1,ndim_kk(me1)
              ir=n_mpi_kk(me1)+ind_k
              t=0.0d0
              t1=0.0d0
              call zgemm('c','c',num_orb_proj,n1,n1,(1.d0,0.d0),
     &          zt(jnd,1,ispin,ind_k),nfun,fifim(1,1,jj,jsort,ispin),
     $          maxel,(0.d0,0.d0),t,num_orb_proj) ! t0=> transpose: usual product basis, and then complex conjugate of product basis
              call zgemm('n','n',num_orb_proj,num_orb_proj,n1,
     $          (1.d0,0.d0),t,num_orb_proj,zt(jnd,1,ispin,ind_k),nfun,
     &          (0.d0,0.d0),t1,num_orb_proj)
              do iorb=1, num_orb_proj ! tau1
                do jorb=1, num_orb_proj ! tau2
                  t2(ir,iorb,jorb)=t1(iorb,jorb) ! 
                enddo
              enddo

              do iat=1, natom_ciw            
                iatom=iatom_ciw(iat)
                isort=is(iatom)
                n=n_ciw_atom(iatom)
                nn=n**2
                do jorb=1, n
                  do iorb=1, n
                    call indexing(nn*nspin,3,(/n,n,nspin/),1,korb,
     $                (/iorb,jorb,ispin/))
                    br_mt(jpb,korb,ind_k,iat)
     $                =t1(ciw_ind(iorb,iat),ciw_ind(jorb,iat)) ! B_{1,2,I}^tau, 1:iorb,2:jorb,:jpb
                  enddo
                enddo
              enddo              
            enddo
            
            call mpi_allreduce_dcmplx(t2,nqdiv*num_orb_proj**2,
     $        mpi_sum,mpi_comm_world)

            call fft3(ndiv(1),ndiv(2),ndiv(3),
     $        num_orb_proj**2,t2,-1)      
            do ind_k=1,ndim3_k(me3_k+1) 
              k=n3_mpi_k(me3_k+1)+ind_k 
              do kst=1,k_star(k)
                iq=k_list(kst,k)
                call zone1_number(pnt(:,iq),rb0,ndiv,ind)   
                do iorb=1, num_orb_proj
                  do jorb=1, num_orb_proj
                    prod_wan2(jpb,iorb,jorb,ispin,ind_k,kst)
     $                =t2(ind,iorb,jorb)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      deallocate(t)
      deallocate(t1)
      deallocate(t2)
      
      
      deallocate(ndim_nr)
      deallocate(n_mpi_nr)
      deallocate(fifim)
      deallocate(pb_int)
      deallocate(zt)
      deallocate(zi)
      deallocate(b0)
      deallocate(b1)
      deallocate(b2)      


      allocate(nwpb(natom_ciw))
      nwpb=0

      do ind_k=1,ndim_kk(me1)
        iq=n_mpi_kk(me1)+ind_k
        k0=i_kref(iq)                                  
        call s_overlap(0,s(1,1,ind_k),nplw_gw,nplwgw(k0),indpw_gw(1,iq))
      enddo

      allocate(pi_ww_c(mx2,mx2*nspin,nspin,natom_ciw))
      pi_ww_c=0.0d0

      allocate(m_p(n_pbtot,mx2*nspin,natom_ciw,ndim_kk(me+1)))
      m_p=0.0d0

!     write(100,*)'n',n,'nplw_gw',nplw_gw,'nplwgw',nplwgw,
!     $ 'br_int',br_int
!     call mpi_bcast(nproc_k,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)

      do iat=1, natom_ciw
        iatom=iatom_ciw(iat)
        isort=is(iatom)
        n=n_ciw_atom(iatom)
        nn=n**2

        allocate(cv_mmmm(nn*nspin,nn*nspin)) ! 
        cv_mmmm=0.0d0
        allocate(cv_mmmm2(nn*nspin,nn*nspin))
        cv_mmmm2=0.0d0        

        do ind_k=1,ndim_kk(me1)
          iq=n_mpi_kk(me1)+ind_k
          k0=i_kref(iq)                            
!     MT region
          call zgemm('c','n',nn*nspin,nn*nspin,
     $      n_pbmt,(1.d0,0.d0),
     &      br_mt(1,1,ind_k,iat),n_pbmt,
     &      br_mt(1,1,ind_k,iat),n_pbmt,        
     &      (1.d0,0.d0),cv_mmmm,nn*nspin) !1,2^*,3^*,4

!     IT region
          
          allocate(tmp(nn*nspin,nplw_gw))
          
          call zgemm('c','n',nn*nspin,nplwgw(k0),
     $      nplwgw(k0),fac0,
     &      br_int(1,1,ind_k,iat),nplw_gw,
     &      s(1,1,ind_k),nplw_gw,(0.d0,0.d0),tmp,nn*nspin)

          call zgemm('n','n',nn*nspin,nn*nspin,
     $      nplwgw(k0),(1.d0,0.d0),
     &      tmp,nn*nspin,br_int(1,1,ind_k,iat),nplw_gw,        
     &      (1.d0,0.d0),cv_mmmm,nn*nspin)
          deallocate(tmp)
        enddo
        if(nproc/=1) call mpi_allreduce_dcmplx
     $    (cv_mmmm,(nn*nspin)**2,mpi_sum,MPI_COMM_WORLD)
        
        ndim=nn*nspin
        
        allocate(diag(ndim))    !
        
        allocate(xx_coef(nn,nn*nspin,nspin))
        
        cv_mmmm2=cv_mmmm
!     write(100,*)'ndim',ndim,'diag',diag,
!     $   'cv_mmmm2',cv_mmmm2 
!     call mpi_bcast(nproc_k,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        
        call hermitianeigen_cmplxdouble(ndim,diag,cv_mmmm2)
        dm=diag(ndim)
        if(maswrk) write(iun,*)
     &    ' WAN_PB_RS: overlap matrix diagonalization'
c     ---- Do not exceed max value of product functions ----------
        j=0
        do i=1,ndim
          if(diag(i)/dm>=eps_pb_vrt) then
            j=j+1
            xx_coef(1:nn,j,1)=cv_mmmm2(1:nn,i)/dsqrt(diag(i))
            if(nspin==2)
     $        xx_coef(1:nn,j,2)=cv_mmmm2((nn+1):ndim,i)/dsqrt(diag(i))
          endif
        enddo
        if(maswrk) write(iun,'(a, i5)')
     $    'iatom=', iatom_ciw(iat)
        if(maswrk) write(iun,'(a,10000000(1x,e12.5))')
     &    ' NPB ', (diag(ndim-i+1), i=1, ndim)       
        if(maswrk) write(iun,'(a,2(i3,2x),2(1x,e12.5))')
     &    ' NPB EIG_max EIG_min  ',j,ndim,diag(ndim),
     &    diag(ndim-j+1)
        nwpb(iat)=j
        deallocate(diag)
c     ------- overlap between product FPB and WPB----

        allocate(tmp(nn,nn))
        do isp=1,nspin
          do jsp=1,nspin
            tmp=0.0d0
            do iorb=1, nn
              do jorb=1, nn
                tmp(iorb,jorb)
     $            =cv_mmmm(iorb+nn*(isp-1),jorb+nn*(jsp-1))
              enddo
            enddo
            call zgemm('n','n',nn,nwpb(iat),nn,(1.d0,0.d0),
     &        tmp,nn,xx_coef(1,1,jsp),nn,
     &        (1.d0,0.d0),pi_ww_c(1,1,isp,iat),mx2) ! (4, 1; wpb), 1 with complex conjugate
          enddo
        enddo
        deallocate(tmp)

!     <M_I^q|Pi_i^q>
!     MT contribution
        
        allocate(tm(n_pbmt,nwpb(iat),ndim_kk(me1))) !
        tm=0.0d0
        do ind_k=1,ndim_kk(me1)
          do ispin=1,nspin
            if (ispin .eq.1) then
              ind=1
            else
              ind=nn+1
            endif
            call zgemm('n','n',n_pbmt,nwpb(iat),nn,(1.d0,0.d0),
     &        br_mt(1,ind,ind_k,iat),n_pbmt,
     &        xx_coef(1,1,ispin),nn,(1.d0,0.d0),tm(1,1,ind_k),
     &        n_pbmt)
          enddo
        enddo
        allocate(tm1(nqdiv,n_pbmt)) !
        do j=1,nwpb(iat)
          tm1=(0.d0,0.d0)
          do i=1,n_pbmt
            do ind_k=1,ndim_kk(me1)
              ir=n_mpi_kk(me1)+ind_k
              tm1(ir,i)=tm(i,j,ind_k)
            enddo
          enddo
          if(nproc/=1) call mpi_allreduce_dcmplx
     $      (tm1,nqdiv*n_pbmt,mpi_sum,MPI_COMM_WORLD)
          call fft3(ndiv(1),ndiv(2),ndiv(3),n_pbmt,tm1,-1)
          do ind_k=1,ndim_kk(me1)
            iq=n_mpi_kk(me1)+ind_k
            call zone1_number(pnt(:,iq),rb0,ndiv,ind)
            do i=1,n_pbmt
              m_p(i,j,iat,ind_k)=tm1(ind,i) ! ind_k zeros bz
            enddo
          enddo
        enddo
        deallocate(tm,tm1)
        
!     interstitial contribution        
        allocate(tm1(nplw_gw,nwpb(iat)))
        do ind_k=1,ndim_kk(me1)
          iq=n_mpi_kk(me1)+ind_k
          iq0=i_kref(iq)
          tm1=(0.d0,0.d0)
          do ispin=1,nspin
            if (ispin .eq.1) then
              ind=1
            else
              ind=nn+1
            endif            
            call zgemm('n','n',nplwgw(iq0),nwpb(iat),nn,(1.d0,0.d0),
     &        br_int(1,ind,ind_k,iat),nplw_gw,
     &        xx_coef(1,1,ispin),nn,(1.d0,0.d0),tm1,nplw_gw)
          enddo

          do i=1,nplwgw(iq0)          
            do j=1,nwpb(iat)
              m_p(i+n_pbmt,j,iat,ind_k)=tm1(i,j)
            enddo
          enddo
        enddo                   !! over ind_k

        deallocate(tm1)
        deallocate(cv_mmmm)
        deallocate(cv_mmmm2)
        deallocate(xx_coef)
      enddo
      
      deallocate(br_mt,br_int)
      deallocate(s)
      
      nwpb_max=maxval(nwpb)
      nwpb_all=sum(nwpb)

      allocate(pi_ww(n_ciw_max,n_ciw_max,nwpb_max,nspin,natom_ciw))
      pi_ww=0.0d0
      allocate(m_pi(n_pbtot,nwpb_max,natom,
     $  ndim3_k(me3_k+1),k_star_max))
      m_pi=0.0d0      
      do iat=1, natom_ciw
        iatom=iatom_ciw(iat)
        n=n_ciw_atom(iatom)
        nn=n**2
        do ispin=1, nspin
          do iwan=1,n           ! 4, without c.c
            do jwan=1,n         ! 1
              call indexing(nn,2,(/n,n/),1,korb,
     $          (/iwan,jwan/))              
              do ipb=1,nwpb(iat) ! without c.c
                pi_ww(iwan,jwan,ipb,ispin,iat)
     $            =pi_ww_c(korb,ipb,ispin,iat)
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(pi_ww_c)
      
      do iat=1, natom_ciw
        iatom=iatom_ciw(iat)
        do ipb=1,n_pbtot
          allocate(tmp(nwpb(iat),nqdiv))
          tmp=0.0d0
          do jpb=1,nwpb(iat)
            do ind_k=1,ndim_kk(me1)
              ik=n_mpi_kk(me1)+ind_k
              tmp(jpb,ik)=m_p(ipb,jpb,iat,ind_k)
            enddo
          enddo
          call mpi_allreduce_dcmplx
     $      (tmp,nwpb(iat)*nqdiv,
     $      mpi_sum,MPI_COMM_WORLD)
          do jpb=1,nwpb(iat)
            do ind_k=1,ndim3_k(me3_k+1) 
              k=n3_mpi_k(me3_k+1)+ind_k 
              do kst=1,k_star(k)
                iq=k_list(kst,k)
                m_pi(ipb,jpb,iatom,ind_k,kst)=tmp(jpb,iq)
              enddo
            enddo
          enddo
          deallocate(tmp)
        enddo
      enddo
      deallocate(m_p)




c$$$  c     -------- MTILD is formed here first ----------------------------

      if (we_need_u) then
c$$$  if (maswrk) then
c$$$  do iat=1, natom_ciw
c$$$  jat=iatom_ciw(iat)
c$$$  write(ch_atom, *) jat
c$$$  allocate(pi_ww_part
c$$$  $        (n_ciw_atom(jat),n_ciw_atom(jat),nwpb(iat),nspin))
c$$$  pi_ww_part=0.0d0
c$$$  pi_ww_part=pi_ww(1:n_ciw_atom(jat),
c$$$  $        1:n_ciw_atom(jat),1:nwpb(iat),:,iat)
c$$$  open(unit=278, file='wannier_prod_bas_'
c$$$  $        //trim(adjustl(ch_atom))//'.dat',
c$$$  $        form = 'unformatted')
c$$$  write(278) n_ciw_atom(jat),nwpb(iat),nspin
c$$$  write(278) pi_ww_part
c$$$  close(278)
c$$$  deallocate(pi_ww_part)
c$$$  enddo
c$$$  endif
        
c$$$  if (me3_tau .eq. 0)then
c$$$  write(179,*) 'prod_wan2_test'
c$$$  allocate(tmp2(num_orb_proj))
c$$$  
c$$$  do ispin=1, nspin
c$$$  do ind_k=1,ndim3_k(me3_k+1) 
c$$$  k=n3_mpi_k(me3_k+1)+ind_k 
c$$$  do kst=1,k_star(k)
c$$$  iq=k_list(kst,k)
c$$$  do ib=1, n_pbtot
c$$$  
c$$$  tmp2=0.0d0
c$$$  
c$$$  do iat=1, natom_ciw
c$$$  iatom=iatom_ciw(iat)
c$$$  do jb=1, nwpb(iat)                    
c$$$  do jorb=1, n_ciw_atom(iatom)
c$$$  iorb=ciw_ind(jorb,iat)
c$$$  
c$$$  tmp2(iorb)=tmp2(iorb)+m_pi(ib,jb,iat,ind_k,kst)
c$$$  $                    *dconjg(pi_ww(jorb,jorb,jb,ispin,iat))
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  
c$$$  write(179,'(4i5, 1000(4f12.6))') ispin,ind_k,kst,ib,
c$$$  $              ((prod_wan2(ib,iorb,iorb,ispin,ind_k,kst),
c$$$  $              tmp2(iorb)), iorb=1, num_orb_proj)
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  deallocate(tmp2)
c$$$  endif
      endif

      pol_iter_for_u=1

      if ((we_need_u) .and. (
     $  (u_mode .eq. 'bnse') .or. (u_mode .eq. 'wnds')      
     $  )) then

c$$$  do ii=1, pol_iter_for_u
c$$$  p_low_string(
        
        if (p_low_mode .eq. 0) then !d_d+p_d+d_p+p_p
          pol_iter_for_u=1
        elseif (p_low_mode .eq. 1) then !d_d
          pol_iter_for_u=1
        elseif (p_low_mode .eq. 2) then ! p_p
          pol_iter_for_u=1          
        elseif (p_low_mode .eq. 3) then !p_d +d_p
          pol_iter_for_u=2                            
        elseif (p_low_mode .eq. 4) then ! d_rh+rh_d
          pol_iter_for_u=2                                          
        elseif (p_low_mode .eq. 5) then ! d_rl+rl_d
          pol_iter_for_u=2                                                        
        elseif (p_low_mode .eq. 6) then !p_rh+rh_p
          pol_iter_for_u=2                                                           
        elseif (p_low_mode .eq. 7) then !p_rl+rl_p
          pol_iter_for_u=2                                                           
        elseif (p_low_mode .eq. 8) then !rh_rl+rl_rh
          pol_iter_for_u=2          
        elseif (p_low_mode .eq. 9) then !all
          pol_iter_for_u=1

          
        elseif (p_low_mode .eq. 12) then ! d_d+p_p
          pol_iter_for_u=2          
        elseif (p_low_mode .eq. 13) then !d_d+p_d +d_p
          pol_iter_for_u=3                            
        elseif (p_low_mode .eq. 14) then ! d_d+d_rh+rh_d
          pol_iter_for_u=3                                          
        elseif (p_low_mode .eq. 15) then !d_d+d_rl+rl_d
          pol_iter_for_u=3                                                        
        elseif (p_low_mode .eq. 16) then !d_d+p_rh+rh_p
          pol_iter_for_u=3                                                           
        elseif (p_low_mode .eq. 17) then !d_d+p_rl+rl_p
          pol_iter_for_u=3                                                           
        elseif (p_low_mode .eq. 18) then !d_d+rh_rl+rl_rh
          pol_iter_for_u=3

        elseif (p_low_mode .eq. 24) then ! d_d+p_d+d_p+p_p+d_rh+rh_d
          pol_iter_for_u=3                                          
        elseif (p_low_mode .eq. 25) then !d_d+p_d+d_p+p_p+d_rl+rl_d
          pol_iter_for_u=3                                                        
        elseif (p_low_mode .eq. 26) then !d_d+p_d+d_p+p_p+p_rh+rh_p
          pol_iter_for_u=3                                                           
        elseif (p_low_mode .eq. 27) then !d_d+p_d+d_p+p_p+p_rl+rl_p
          pol_iter_for_u=3                                                           
        elseif (p_low_mode .eq. 28) then !d_d+p_d+d_p+p_p+rh_rl+rl_rh
          pol_iter_for_u=3

        elseif (p_low_mode .eq. 35) then !d_d+p_d+d_p+p_p+d_rl+rl_d
          pol_iter_for_u=4                                                        
        elseif (p_low_mode .eq. 36) then !d_d+p_d+d_p+p_p+p_rh+rh_p
          pol_iter_for_u=4                                                           
        elseif (p_low_mode .eq. 37) then !d_d+p_d+d_p+p_p+p_rl+rl_p
          pol_iter_for_u=4                                                           
        elseif (p_low_mode .eq. 38) then !d_d+p_d+d_p+p_p+rh_rl+rl_rh
          pol_iter_for_u=4                                                                               
c$$$  elseif (p_low_mode .eq. 9) then !all
c$$$  pol_iter_for_u=1                                                           
        endif

        allocate(include_bands_foru(nbndf,ndim3_k(me3_k+1),ispin,2,
     $    pol_iter_for_u))      !it,prtausolid_index
        include_bands_foru=0
        
        do ind_k=1,ndim3_k(me3_k+1)
          ik=n3_mpi_k(me3_k+1)+ind_k
          do ispin=1,nspin
            
            if (p_low_mode .eq. 0) then !d_d +p_d +d_p+p_p
              do ibnd=1,n_bnd(ik,ispin)              
                if ((ind_bands_foru(ibnd,ik,ispin) .eq. 1)
     $            .or. (ind_bands_foru(ibnd,ik,ispin) .eq. 2)) then
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif
              enddo
              
            elseif (p_low_mode .eq. 1) then !d_d
              do ibnd=1,n_bnd(ik,ispin)              
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 1) then
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif
              enddo

            elseif (p_low_mode .eq.2) then ! p_p
              do ibnd=1,n_bnd(ik,ispin)
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 2) then
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif
              enddo
              
            elseif (p_low_mode .eq.3) then !p_d +d_p
              do ibnd=1,n_bnd(ik,ispin)

                
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 1) then                
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,2)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 2) then
                  include_bands_foru(ibnd,ind_k,ispin,1,2)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif
              enddo

            elseif (p_low_mode .eq. 4) then ! d_rh+rh_d
              do ibnd=1,n_bnd(ik,ispin)

                if (ind_bands_foru(ibnd,ik,ispin) .eq. 1) then                
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,2)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 3) then
                  include_bands_foru(ibnd,ind_k,ispin,1,2)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif
              enddo


            elseif (p_low_mode .eq.5) then !d_rl+rl_d
              do ibnd=1,n_bnd(ik,ispin)

                if (ind_bands_foru(ibnd,ik,ispin) .eq. 1) then                
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,2)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 4) then
                  include_bands_foru(ibnd,ind_k,ispin,1,2)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif
              enddo



            elseif (p_low_mode .eq.6) then !p_rh+rh_p

              do ibnd=1,n_bnd(ik,ispin)

                if (ind_bands_foru(ibnd,ik,ispin) .eq. 2) then                
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,2)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 3) then
                  include_bands_foru(ibnd,ind_k,ispin,1,2)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif
              enddo


            elseif (p_low_mode .eq.7) then !p_rl+rl_p

              do ibnd=1,n_bnd(ik,ispin)

                if (ind_bands_foru(ibnd,ik,ispin) .eq. 2) then                
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,2)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 4) then
                  include_bands_foru(ibnd,ind_k,ispin,1,2)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif
              enddo

            elseif (p_low_mode .eq.8) then !rh_rl+rl_rh
              do ibnd=1,n_bnd(ik,ispin)
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 3) then
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,2)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 4) then
                  include_bands_foru(ibnd,ind_k,ispin,1,2)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif
              enddo              

            elseif (p_low_mode .eq. 9) then !all
              do ibnd=1,n_bnd(ik,ispin)
                include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                include_bands_foru(ibnd,ind_k,ispin,2,1)=1
              enddo
              
            elseif (p_low_mode .eq.12) then ! d_d+p_p
              do ibnd=1,n_bnd(ik,ispin)
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 1) then
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif                
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 2) then
                  include_bands_foru(ibnd,ind_k,ispin,1,2)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,2)=1
                endif
              enddo
              
            elseif (p_low_mode .eq.13) then !d_d+p_d +d_p
              do ibnd=1,n_bnd(ik,ispin)
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 1) then
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif
                
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 1) then                
                  include_bands_foru(ibnd,ind_k,ispin,1,2)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,3)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 2) then
                  include_bands_foru(ibnd,ind_k,ispin,1,3)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,2)=1
                endif
              enddo

            elseif (p_low_mode .eq. 14) then ! d_d+d_rh+rh_d
              do ibnd=1,n_bnd(ik,ispin)
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 1) then
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 1) then                
                  include_bands_foru(ibnd,ind_k,ispin,1,2)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,3)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 3) then
                  include_bands_foru(ibnd,ind_k,ispin,1,3)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,2)=1
                endif
              enddo


            elseif (p_low_mode .eq.15) then !d_d+d_rl+rl_d
              do ibnd=1,n_bnd(ik,ispin)
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 1) then
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 1) then                
                  include_bands_foru(ibnd,ind_k,ispin,1,2)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,3)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 4) then
                  include_bands_foru(ibnd,ind_k,ispin,1,3)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,2)=1
                endif
              enddo



            elseif (p_low_mode .eq.16) then !d_d+p_rh+rh_p

              do ibnd=1,n_bnd(ik,ispin)
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 1) then
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif                
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 2) then                
                  include_bands_foru(ibnd,ind_k,ispin,1,2)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,3)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 3) then
                  include_bands_foru(ibnd,ind_k,ispin,1,3)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,2)=1
                endif
              enddo


            elseif (p_low_mode .eq.17) then !d_d+p_rl+rl_p

              do ibnd=1,n_bnd(ik,ispin)
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 1) then
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif                
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 2) then                
                  include_bands_foru(ibnd,ind_k,ispin,1,2)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,3)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 4) then
                  include_bands_foru(ibnd,ind_k,ispin,1,3)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,2)=1
                endif
              enddo

            elseif (p_low_mode .eq.18) then !d_d+rh_rl+rl_rh
              do ibnd=1,n_bnd(ik,ispin)
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 1) then
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif                
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 3) then
                  include_bands_foru(ibnd,ind_k,ispin,1,2)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,3)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 4) then
                  include_bands_foru(ibnd,ind_k,ispin,1,3)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,2)=1
                endif
              enddo




            elseif (p_low_mode .eq. 24) then ! d_d+d_rh+rh_d
              do ibnd=1,n_bnd(ik,ispin)
                if ((ind_bands_foru(ibnd,ik,ispin) .eq. 1)
     $            .or. (ind_bands_foru(ibnd,ik,ispin) .eq. 2)) then
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 1) then                
                  include_bands_foru(ibnd,ind_k,ispin,1,2)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,3)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 3) then
                  include_bands_foru(ibnd,ind_k,ispin,1,3)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,2)=1
                endif
              enddo


            elseif (p_low_mode .eq.25) then !d_d+d_rl+rl_d
              do ibnd=1,n_bnd(ik,ispin)
                if ((ind_bands_foru(ibnd,ik,ispin) .eq. 1)
     $            .or. (ind_bands_foru(ibnd,ik,ispin) .eq. 2)) then
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 1) then                
                  include_bands_foru(ibnd,ind_k,ispin,1,2)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,3)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 4) then
                  include_bands_foru(ibnd,ind_k,ispin,1,3)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,2)=1
                endif
              enddo

            elseif (p_low_mode .eq.26) then !d_d+p_rh+rh_p

              do ibnd=1,n_bnd(ik,ispin)
                if ((ind_bands_foru(ibnd,ik,ispin) .eq. 1)
     $            .or. (ind_bands_foru(ibnd,ik,ispin) .eq. 2)) then
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 2) then                
                  include_bands_foru(ibnd,ind_k,ispin,1,2)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,3)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 3) then
                  include_bands_foru(ibnd,ind_k,ispin,1,3)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,2)=1
                endif
              enddo


            elseif (p_low_mode .eq.27) then !d_d+p_rl+rl_p

              do ibnd=1,n_bnd(ik,ispin)
                if ((ind_bands_foru(ibnd,ik,ispin) .eq. 1)
     $            .or. (ind_bands_foru(ibnd,ik,ispin) .eq. 2)) then
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif

                if (ind_bands_foru(ibnd,ik,ispin) .eq. 2) then                
                  include_bands_foru(ibnd,ind_k,ispin,1,2)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,3)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 4) then
                  include_bands_foru(ibnd,ind_k,ispin,1,3)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,2)=1
                endif
              enddo

            elseif (p_low_mode .eq.28) then !d_d+rh_rl+rl_rh
              do ibnd=1,n_bnd(ik,ispin)
                if ((ind_bands_foru(ibnd,ik,ispin) .eq. 1)
     $            .or. (ind_bands_foru(ibnd,ik,ispin) .eq. 2)) then
                  include_bands_foru(ibnd,ind_k,ispin,1,1)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,1)=1
                endif

                if (ind_bands_foru(ibnd,ik,ispin) .eq. 3) then
                  include_bands_foru(ibnd,ind_k,ispin,1,2)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,3)=1
                endif
                if (ind_bands_foru(ibnd,ik,ispin) .eq. 4) then
                  include_bands_foru(ibnd,ind_k,ispin,1,3)=1
                  include_bands_foru(ibnd,ind_k,ispin,2,2)=1
                endif
              enddo              

            endif
          enddo
        enddo
      elseif ((we_need_u) .and. (
     $    (u_mode .eq. 'edmf')
     $    )) then
        pol_iter_for_u=0
      endif

      if ((we_need_u) .and. (u_mode .eq. 'wnds')) then
        deallocate(eigenvalues_cor)
        deallocate(v_matrix_cor)
        if (num_orb_iter .gt. 0) then
          deallocate(eigenvalues_iter)
          deallocate(v_matrix_iter)
        endif

        deallocate(eigenvalues_rest)
        deallocate(v_matrix_rest)

        deallocate(sorting_disen)
        deallocate(eval_disen_sort)
      endif      

      


      allocate(w_k_den(
     $  num_orb_proj,nspin,num_orb_proj,nspin,k_star_max,
     $  ndim3_k(me3_k+1), ndim3_nu))
      w_k_den=0.0d0
      allocate(w_loc(
     $  num_orb_proj,num_orb_proj,num_orb_proj,num_orb_proj,nspin,nspin,
     $  ndim3_nu))
      w_loc=0.0d0

      allocate(wc_k_den_tau(
     $  num_orb_proj,nspin,num_orb_proj,nspin,k_star_max,
     $  ndim3_k(me3_k+1), ndim3_tau))
      wc_k_den_tau=0.0d0            

      allocate(epsinv_k_den(
     $  num_orb_proj,nspin,num_orb_proj,nspin,k_star_max,
     $  ndim3_k(me3_k+1), ndim3_nu))
      epsinv_k_den=0.0d0
      allocate(epsinv_loc(
     $  num_orb_proj,num_orb_proj,num_orb_proj,num_orb_proj,nspin,nspin,
     $  ndim3_nu))
      epsinv_loc=0.0d0

      allocate(eps_k_den(
     $  num_orb_proj,nspin,num_orb_proj,nspin,k_star_max,
     $  ndim3_k(me3_k+1), ndim3_nu))
      eps_k_den=0.0d0
      allocate(eps_loc(
     $  num_orb_proj,num_orb_proj,num_orb_proj,num_orb_proj,nspin,nspin,
     $  ndim3_nu))
      eps_loc=0.0d0

      allocate(chi0_k_den(
     $  num_orb_proj,nspin,num_orb_proj,nspin,k_star_max,
     $  ndim3_k(me3_k+1), ndim3_nu))
      chi0_k_den=0.0d0
      allocate(chi0_loc(
     $  num_orb_proj,num_orb_proj,num_orb_proj,num_orb_proj,nspin,nspin,
     $  ndim3_nu))
      chi0_loc=0.0d0

      allocate(chi_k_den(
     $  num_orb_proj,nspin,num_orb_proj,nspin,k_star_max,
     $  ndim3_k(me3_k+1), ndim3_nu))
      chi_k_den=0.0d0
      allocate(chi_k_den_tau(
     $  num_orb_proj,nspin,num_orb_proj,nspin,k_star_max,
     $  ndim3_k(me3_k+1), ndim3_tau))
      chi_k_den_tau=0.0d0

      allocate(chi_loc(
     $  num_orb_proj,num_orb_proj,num_orb_proj,num_orb_proj,nspin,nspin,
     $  ndim3_nu))
      chi_loc=0.0d0

      allocate(u_k_den(
     $  num_orb_proj,nspin,num_orb_proj,nspin,k_star_max,
     $  ndim3_k(me3_k+1), ndim3_nu))
      u_k_den=0.0d0

      allocate(v_k_den(
     $  num_orb_proj,nspin,num_orb_proj,nspin,k_star_max,
     $  ndim3_k(me3_k+1)))
      v_k_den=0.0d0      
      

      allocate(w_wpb_k(
     $  nwpb_max,natom_ciw,nwpb_max,natom_ciw,k_star_max,
     $  ndim3_k(me3_k+1), ndim3_nu))
      w_wpb_k=0.0d0

      allocate(chi0_wpb_k(
     $  nwpb_max,natom_ciw,nwpb_max,natom_ciw,k_star_max,
     $  ndim3_k(me3_k+1), ndim3_nu))
      chi0_wpb_k=0.0d0

      allocate(w_wpb_local(
     $  nwpb_max,natom_ciw,nwpb_max,natom_ciw,ndim3_nu))
      w_wpb_local=0.0d0

      allocate(chi0_wpb_local(
     $  nwpb_max,natom_ciw,nwpb_max,natom_ciw,ndim3_nu))
      chi0_wpb_local=0.0d0
      
      allocate(v_loc(
     $  num_orb_proj,num_orb_proj,num_orb_proj,num_orb_proj,
     $  nspin,nspin))
      v_loc=0.0d0        
      

      if (pol_iter_for_u .gt. 0) then
        
        allocate(u_loc(
     $    num_orb_proj,num_orb_proj,num_orb_proj,num_orb_proj,
     $    nspin,nspin,ndim3_nu))
        u_loc=0.0d0

        
        
        allocate(chi0_low_k_den(
     $    num_orb_proj,nspin,num_orb_proj,nspin,k_star_max,
     $    ndim3_k(me3_k+1), ndim3_nu))
        chi0_low_k_den=0.0d0
        allocate(chi0_low_loc(
     $    num_orb_proj,num_orb_proj,num_orb_proj,num_orb_proj,
     $    nspin,nspin,ndim3_nu))
        chi0_low_loc=0.0d0
        
        allocate(epsinv_high_k_den(
     $    num_orb_proj,nspin,num_orb_proj,nspin,k_star_max,
     $    ndim3_k(me3_k+1), ndim3_nu))
        epsinv_high_k_den=0.0d0
        allocate(epsinv_high_loc(
     $    num_orb_proj,num_orb_proj,num_orb_proj,num_orb_proj,
     $    nspin,nspin,ndim3_nu))
        epsinv_high_loc=0.0d0
        
        allocate(eps_high_k_den(
     $    num_orb_proj,nspin,num_orb_proj,nspin,k_star_max,
     $    ndim3_k(me3_k+1), ndim3_nu))
        eps_high_k_den=0.0d0
        allocate(eps_high_loc(
     $    num_orb_proj,num_orb_proj,num_orb_proj,num_orb_proj,
     $    nspin,nspin,ndim3_nu))
        eps_high_loc=0.0d0
        
        
        allocate(chi0_high_k_den(
     $    num_orb_proj,nspin,num_orb_proj,nspin,k_star_max,
     $    ndim3_k(me3_k+1), ndim3_nu))
        chi0_high_k_den=0.0d0
        allocate(chi0_high_loc(
     $    num_orb_proj,num_orb_proj,num_orb_proj,num_orb_proj,
     $    nspin,nspin,ndim3_nu))
        chi0_high_loc=0.0d0
        
        allocate(chi_high_k_den(
     $    num_orb_proj,nspin,num_orb_proj,nspin,k_star_max,
     $    ndim3_k(me3_k+1), ndim3_nu))
        chi_high_k_den=0.0d0
        
        allocate(chi_high_k_den_tau(
     $    num_orb_proj,nspin,num_orb_proj,nspin,k_star_max,
     $    ndim3_k(me3_k+1), ndim3_tau))
        chi_high_k_den_tau=0.0d0
        
        allocate(chi_high_loc(
     $    num_orb_proj,num_orb_proj,num_orb_proj,num_orb_proj,
     $    nspin,nspin,ndim3_nu))
        chi_high_loc=0.0d0
      endif
      
c$$$      if (u_mode .eq. 'edmf') then
c$$$        call chi0_projection      
c$$$c$$$        call u_weiss
c$$$      endif

C$$$  allocate(ndim_wpbl(nproc_k))
C$$$  ndim_wpbl=0
C$$$  allocate(n_mpi_wpbl(nproc_k))
C$$$  n_mpi_wpbl=0
c$$$  if (irel .le.1) then
C$$$  call size_shift_par((num_orb_proj*nspin)**2, nproc_k,
C$$$  $  ndim_wpbl,n_mpi_wpbl)
c$$$  else
c$$$  call size_shift_par(((num_orb_cor/2)*nspin)**2,
c$$$  $    nproc_k,ndim_wpbl,n_mpi_wpbl)        
c$$$  endif


C$$$  if (we_need_u) then
C$$$  if (maswrk) then
C$$$  open(unit=278, file='wannier_prod_bas.dat',
C$$$  $      form = 'formatted')
C$$$  write(278,*) mx2, nwp_max, nspin, natom_ciw
C$$$  do iat=1, natom_ciw
C$$$  do kk=1, nspin
C$$$  do jj=1, nwpb(iat)
C$$$  do ii=1, mx2                
C$$$  write(278,'(2f20.15)') pi_ww(ii,jj,kk,iat)
C$$$  enddo
C$$$  enddo
C$$$  enddo
C$$$  enddo
C$$$  close(278)
C$$$  endif
C$$$  endif

C$$$  if(allocated(den_wanpb)) deallocate(den_wanpb)
C$$$  allocate(den_wanpb(n_ciw_max,nwp_max,nspin,natom_ciw))
C$$$  den_wanpb=0.0d0

C$$$  !     ! here
      

C$$$  do isp=1,nspin
C$$$  do iat=1, natom_ciw
C$$$  iatom=iatom_ciw(iat)
C$$$  do aa=1,n_ciw_atom(iatom)
C$$$  m12=(aa-1)*n_ciw_atom(iatom)+aa
C$$$  do kk=1, nwpb(iat)
C$$$  den_wanpb(aa,kk,isp,iat)
C$$$  $          =den_wanpb(aa,kk,isp,iat)
C$$$  $          +pi_ww(m12,kk,isp,iat)
C$$$  enddo
C$$$  enddo
C$$$  enddo
C$$$  enddo
c$$$  else
c$$$  
c$$$  
c$$$  do isp=1,nspin
c$$$  do li=1, nfun0_loc(isort)/2
c$$$  l=li_loc(li*2-1,isort)/2
c$$$  mdim=4*l+2
c$$$  
c$$$  allocate(rotmat_mlms2jmj(mdim, mdim))
c$$$  rotmat_mlms2jmj=0.0d0
c$$$  allocate(rotmat_cmplx2real(mdim/2, mdim/2))
c$$$  rotmat_cmplx2real=0.0d0
c$$$  allocate(rotmat_temp(mdim, mdim))
c$$$  rotmat_temp=0.0d0
c$$$  allocate(rotmat_temp2(mdim, mdim))
c$$$  rotmat_temp2=0.0d0                                                
c$$$  call cal_rotmat_mlms2jmj(l, rotmat_mlms2jmj)
c$$$  call cal_rotmat_cmplx2real
c$$$  $          (mdim/2,rotmat_cmplx2real)
c$$$  rotmat_temp(1:(mdim/2), 1:(mdim/2))
c$$$  $          =rotmat_cmplx2real
c$$$  rotmat_temp((mdim/2+1):(mdim),
c$$$  $          (mdim/2+1):(mdim))
c$$$  $          =rotmat_cmplx2real
c$$$  rotmat_temp2=matmul(rotmat_temp,
c$$$  $          dconjg(transpose(rotmat_mlms2jmj)))
c$$$  do aa=1,mdim/2
c$$$  do ispin=1, 2
c$$$  aaindex=aa+(ispin-1)*(mdim/2)
c$$$  
c$$$  ind1=0
c$$$  do li_i=2*li-1, 2*li
c$$$  do lim_i=1,(li_loc(li_i,isort)+1)/2*2
c$$$  ind1=ind1+1
c$$$  m1=indbasa_loc(lim_i,li_i,isort)
c$$$  
c$$$  ind2=0
c$$$  do li_j=2*li-1, 2*li
c$$$  do lim_j=1, (li_loc(li_j,isort)+1)/2*2
c$$$  ind2=ind2+1
c$$$  m2=indbasa_loc(lim_j,li_j,isort)
c$$$  
c$$$  m12=(m2-1)*nfun_loc(isort)+m1                      
c$$$  do kk=1, nwpb(iatom_c)
c$$$  den_wanpb(li,isp,kk,iatom_c)
c$$$  $                        =den_wanpb(li,isp,kk,iatom_c)
c$$$  $                        +rotmat_temp2(aaindex,ind1)
c$$$  $                        *dconjg(rotmat_temp2(aaindex,ind2))
c$$$  $                        *pi_ww(m12,kk,isp,iatom_c)/dble(mdim)
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  deallocate(rotmat_mlms2jmj)
c$$$  deallocate(rotmat_cmplx2real)
c$$$  deallocate(rotmat_temp)
c$$$  deallocate(rotmat_temp2)
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  endif 


c$$$  
c$$$  do iat=1,natom_ciw
c$$$  iatom=iatom_ciw(iat)
c$$$  isort=is(iatom)
c$$$  
c$$$  allocate(ttt1(n_pbtot,nwpb(iat),k_star_max))          
c$$$  allocate(ttt2(n_pbtot,nwpb(iat),k_star_max))
c$$$  allocate(tmp(nwpb,nwpb))        
c$$$  
c$$$  do ind_k=1,ndim3_k(me3_k+1)
c$$$  k=n3_mpi_k(me3_k+1)+ind_k        
c$$$  ttt1=0.0d0
c$$$  ttt2=0.0d0
c$$$  do kst=1,k_star(k)                        
c$$$  do ind_b=1,ndim_pbt(me3_tau+1)
c$$$  i_pb=n_mpi_pbt(me3_tau+1)+ind_b
c$$$  do iw=1, nwpb(iat)
c$$$  ttt1(i_pb,iw,kst)
c$$$  $            =mtild_pi(ind_b,iw,kst,ind_k)            
c$$$  ttt2(i_pb,iw,kst)
c$$$  $            =m_pi(ind_b,iw,kst,ind_k)              
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  
c$$$  if(nproc_tau/=1) then
c$$$  call mpi_allreduce_dcmplx
c$$$  $      (ttt1,n_pbtot*nwpb*k_star_max,mpi_sum,comm_k)
c$$$  call mpi_allreduce_dcmplx
c$$$  $      (ttt2,n_pbtot*nwpb*k_star_max,mpi_sum,comm_k)
c$$$  endif      
c$$$  
c$$$  
c$$$  
c$$$  do kst=1,k_star(k)
c$$$  iq=k_list(kst,k)
c$$$  tmp=0.0d0
c$$$  
c$$$  call zgemm('c','n',nwpb,nwpb,
c$$$  $      n_pbtot,(1.d0,0.d0), ttt1(1,1,kst),n_pbtot,
c$$$  &      ttt2(1,1,kst),n_pbtot,(0.d0,0.d0),tmp,nwpb)
c$$$  write(179,*) 'kpoint', ind_k,kst,iq
c$$$  c$$$  call cal_identity_error_dcmplx(tmp,nwpb)
c$$$  
c$$$  c$$$  tmp=0.0d0
c$$$  c$$$  call zgemm('c','n',nwpb,nwpb,
c$$$  c$$$  $      n_pbmt,(1.d0,0.d0), ttt1(1,1,kst),n_pbtot,
c$$$  c$$$  &      ttt2(1,1,kst),n_pbtot,(0.d0,0.d0),tmp,nwpb)
c$$$  c$$$  write(179,*) 'kpoint', ind_k,kst,iq
c$$$  c$$$  call cal_identity_error_dcmplx(tmp,nwpb)
c$$$  enddo
c$$$  enddo
c$$$  deallocate(ttt1,ttt2,tmp)      
      end

      subroutine projector_compare(vmat)

      use parallel_mod
      use solid_mod      
      use wannier90_mod
      
      implicit none
      include "mpif.h"
      complex*16, intent(in) ::
     $  vmat(num_bands,num_orb_proj,ndim_kk(me+1))
      
      integer :: ik,irk,iorb,ibnd,jorb,iprj,jprj,iatom,icell,
     $  jcell,jatom,
     $  indtmp,jndtmp,iband,ii,me1
      double precision :: inv_numk
      complex*16, allocatable :: overlap_mat(:,:),tempmat1_k(:,:),
     $  tempmat2_k(:,:)

      allocate(overlap_mat(num_orb_proj, num_orb_proj))
      overlap_mat=0.0d0
      allocate(tempmat1_k(num_bands,num_orb_proj))
      tempmat1_k=0.0d0
      allocate(tempmat2_k(num_bands,num_orb_proj))
      tempmat2_k=0.0d0      

      inv_numk=1.0d0/dble(nqdiv)

      overlap_mat=0.0d0
      me1=me+1

      do ik=1, ndim_kk(me1)
        irk=n_mpi_kk(me1)+ik
        tempmat1_k=0.0d0
        tempmat2_k=0.0d0        
        do iorb=1, num_orb_cor
          do ibnd=1, num_bands
            tempmat1_k(ibnd,iorb)
     $        =v_matrix(ibnd,wan_correlated(iorb),irk)
            tempmat2_k(ibnd,iorb)=vmat(ibnd,iorb,ik)
          enddo
        enddo

        if (num_orb_iter .gt. 0) then
          do iorb=1, num_orb_iter
            do ibnd=1, num_bands
              tempmat1_k(ibnd,iorb+num_orb_cor)
     $          =v_matrix(ibnd,wan_itinerant(iorb),irk)
              tempmat2_k(ibnd,iorb+num_orb_cor)
     $          =vmat(ibnd,iorb+num_orb_cor,ik)
            enddo
          enddo
        endif

        call zgemm('c','n',num_orb_proj,num_orb_proj,num_bands,
     $    (1.0d0,0.d0),tempmat1_k,num_bands,tempmat2_k,num_bands,
     &    (1.d0,0.d0),overlap_mat,num_orb_proj)            
      enddo      
      
      call mpi_allreduce_dcmplx(overlap_mat,(num_orb_proj)**2,
     $  mpi_sum,mpi_comm_world)
      
      if (maswrk) then
        write(993,*)
     $    "overlap_mat of new projector and old projector"
        do iprj=1, num_orb_proj
          write(993,'(100(2f12.6,2x))')
     $      (overlap_mat(iprj,jprj)*inv_numk, jprj=1, num_orb_proj)
        enddo
      endif

      end


      subroutine projector_renormalization(vmat)

      use parallel_mod
      use wannier90_mod
      use units_mod
      use manager_mod            
      
      implicit none

      complex*16, intent(out) ::
     $  vmat(num_bands,num_orb_proj,ndim_kk(me+1))
      integer :: ik,irk,ibnd,iorb,jorb,iatom,ind_temp,icell,ii,jj,
     $  me1,ndim_kk_max,indtmp,iwan,jbnd


      double precision, allocatable ::
     $  diag_cor(:), diag_iter(:),eig_opt(:)
      complex*16, allocatable :: proj_k_cor(:,:),proj_k_iter(:,:),
     $  proj_k_temp_cor(:,:),proj_k_temp_iter(:,:),
     $  overlapmat_cor(:,:),overlapmat_iter(:,:),
     $  t0_cor(:,:),t0_iter(:,:),
     $  t1_cor(:,:),t1_iter(:,:),tempmat(:,:),tempmat2(:,:)


      me1=me+1
      ndim_kk_max=maxval(ndim_kk)
      
      vmat=0.0d0

      allocate(eig_opt(num_wann))

      allocate(tempmat(num_wann,num_wann))
      tempmat=0.0d0
      allocate(tempmat2(num_wann,num_orb_proj))
      tempmat=0.0d0            
      
      allocate(diag_cor(num_orb_cor)) !
      diag_cor=0.0d0
      allocate(proj_k_cor(num_wann, num_orb_cor)) !
      proj_k_cor=0.0d0
      allocate(proj_k_temp_cor(num_wann,num_orb_cor)) !
      proj_k_temp_cor=0.0d0
      allocate(overlapmat_cor(num_orb_cor,num_orb_cor)) !
      overlapmat_cor=0.0d0
      allocate(t0_cor(num_orb_cor,num_orb_cor)) !
      t0_cor=0.0d0
      allocate(t1_cor(num_orb_cor,num_orb_cor)) !
      t1_cor=0.0d0


      if (num_orb_iter .gt. 0) then
        allocate(diag_iter(num_orb_iter)) !
        diag_iter=0.0d0
        allocate(proj_k_iter(num_wann, num_orb_iter)) !
        proj_k_iter=0.0d0
        allocate(proj_k_temp_iter(num_wann,num_orb_iter)) !
        proj_k_temp_iter=0.0d0
        allocate(overlapmat_iter(num_orb_iter,num_orb_iter)) !
        overlapmat_iter=0.0d0
        allocate(t0_iter(num_orb_iter,num_orb_iter))
        t0_iter=0.0d0
        allocate(t1_iter(num_orb_iter,num_orb_iter))
        t1_iter=0.0d0
      endif

      
      

      do ik=1, ndim_kk(me1)
        irk=n_mpi_kk(me1)+ik
        proj_k_cor=0.0d0
        proj_k_temp_cor=0.0d0
        
        eig_opt=0.0d0
        indtmp=0        
        do ibnd=1, num_bands
          if (lwindow(ibnd,irk)) then
            indtmp=indtmp+1                      
            do iwan=1, num_wann
              eig_opt(iwan)=
     $          cdabs(u_matrix_opt(indtmp,iwan,irk))**2
     $          *eigenvalues(ibnd,irk)
            enddo
          endif
        enddo
        
        
        do ibnd=1, num_wann
          if ((eig_opt(ibnd) .gt. projemin*rydberg) .and.
     $      (eig_opt(ibnd) .lt. projemax*rydberg)) then
            do iorb=1, num_orb_cor
              proj_k_cor(ibnd,iorb)
     $          =proj_k_cor(ibnd,iorb)
     $          +u_matrix(ibnd,wan_correlated(iorb),irk)                    
            enddo
          endif
        enddo
        overlapmat_cor=0.0d0
        diag_cor=0.0d0
        t0_cor=0.0d0
        t1_cor=0.0d0
        call zgemm('c','n',num_orb_cor,num_orb_cor,num_wann,
     $    (1.d0,0.d0),
     $    proj_k_cor,num_wann,
     $    proj_k_cor,num_wann,
     $    (0.d0,0.d0),
     $    overlapmat_cor,num_orb_cor)
C$$$  call eig_val_solver(num_orb_cor,num_orb_cor,overlapmat,diag)
        call hermitianeigen_cmplxdouble
     $    (num_orb_cor,diag_cor,overlapmat_cor)
c     ------ Here we get T = O^(-1/2) -----------------------------------
c        write(179,*) 'in projector_renormalization'
c        write(179,*) 'eigenvalues of overlap matrix'
c        write(179,*) diag_cor
        do iorb=1,num_orb_cor
          t0_cor(:,iorb)=overlapmat_cor(:,iorb)/dsqrt(diag_cor(iorb))
        enddo
        call zgemm('n','c',num_orb_cor,num_orb_cor,num_orb_cor,
     $    (1.d0,0.d0),t0_cor,num_orb_cor,overlapmat_cor,num_orb_cor,
     $    (0.d0,0.d0),t1_cor,num_orb_cor)
        
c$$$  --------- Orthonormalized projector ----------------------
        call zgemm('n','n',num_wann,num_orb_cor,num_orb_cor,
     $    (1.d0,0.d0),proj_k_cor,num_wann,t1_cor,num_orb_cor,
     $    (0.d0,0.d0),proj_k_temp_cor,num_wann)

        tempmat2=0.0d0
        do iorb=1, num_orb_cor
          do ibnd=1, num_wann
            tempmat2(ibnd,iorb)
     $        =tempmat2(ibnd,iorb)
     $        +proj_k_temp_cor(ibnd,iorb)
          enddo
        enddo        


        if (num_orb_iter .gt. 0) then
          proj_k_iter=0.0d0
          proj_k_temp_iter=0.0d0
          
          tempmat=0.0d0
          
          call zgemm('n','c',num_wann,num_wann,num_orb_cor,
     $      (-1.d0,0.d0),proj_k_temp_cor,num_wann,
     &      proj_k_temp_cor,num_wann,(0.d0,0.d0),
     $      tempmat,num_wann)
          do ibnd=1, num_wann
            tempmat(ibnd,ibnd)=1.0d0+tempmat(ibnd,ibnd)
          enddo
          
          do ibnd=1, num_wann
            do jbnd=1, num_wann
              do iorb=1, num_orb_iter
                proj_k_iter(ibnd,iorb)
     $            =proj_k_iter(ibnd,iorb)
     $            +tempmat(ibnd,jbnd)
     $            *u_matrix(jbnd,wan_itinerant(iorb),irk)                    
              enddo
            enddo
          enddo
!     overlap matrix
          overlapmat_iter=0.0d0
          diag_iter=0.0d0
          t0_iter=0.0d0
          t1_iter=0.0d0
          call zgemm('c','n',num_orb_iter,num_orb_iter,num_wann,
     $      (1.d0,0.d0),
     $      proj_k_iter,num_wann,
     $      proj_k_iter,num_wann,
     $      (0.d0,0.d0),
     $      overlapmat_iter,num_orb_iter)
C$$$  call eig_val_solver(num_orb_iter,num_orb_iter,overlapmat,diag)
          call hermitianeigen_cmplxdouble
     $      (num_orb_iter,diag_iter,overlapmat_iter)
c     ------ Here we get T = O^(-1/2) -----------------------------------
c          write(179,*) 'in projector_renormalization'
c          write(179,*) 'eigenvalues of overlap matrix'
c          write(179,*) diag_iter
          do iorb=1,num_orb_iter
            t0_iter(:,iorb)=overlapmat_iter(:,iorb)
     $        /dsqrt(diag_iter(iorb))
          enddo
          call zgemm('n','c',num_orb_iter,num_orb_iter,num_orb_iter,
     $      (1.d0,0.d0),t0_iter,num_orb_iter,overlapmat_iter,
     $      num_orb_iter,(0.d0,0.d0),t1_iter,num_orb_iter)
          
c$$$  --------- Orthonormalized projector ----------------------
          call zgemm('n','n',num_wann,num_orb_iter,num_orb_iter,
     $      (1.d0,0.d0),proj_k_iter,num_wann,t1_iter,num_orb_iter,
     $      (0.d0,0.d0),proj_k_temp_iter,num_wann)
          
          do iorb=1, num_orb_iter
            do ibnd=1, num_wann
              tempmat2(ibnd,iorb+num_orb_cor)
     $          =tempmat2(ibnd,iorb+num_orb_cor)
     $          +proj_k_temp_iter(ibnd,iorb)
            enddo
          enddo
        endif
        
        indtmp=0
        do jbnd=1, num_bands
          if (lwindow(jbnd,irk)) then
            indtmp=indtmp+1            
            do iorb=1, num_orb_proj
              do ibnd=1, num_wann
                vmat(jbnd,iorb,ik)
     $            =vmat(jbnd,iorb,ik)
     $            +u_matrix_opt(indtmp,ibnd,irk)*tempmat2(ibnd,iorb)
              enddo
            enddo
          endif
        enddo        

      enddo

      if (num_orb_iter .gt. 0) then
        deallocate(diag_iter)
        deallocate(proj_k_iter)
        deallocate(proj_k_temp_iter)
        deallocate(overlapmat_iter)
        deallocate(t0_iter)
        deallocate(t1_iter)
      endif
      deallocate(diag_cor)
      deallocate(proj_k_cor)
      deallocate(proj_k_temp_cor)
      deallocate(overlapmat_cor)
      deallocate(t0_cor)
      deallocate(t1_cor)
      deallocate(eig_opt)
      deallocate(tempmat)
      

      end      
      

      subroutine rotate_umatrix(norb,umatrix, rotmat,daggerflag)
      implicit none

      integer, intent(in) :: norb,daggerflag
      complex*16, intent(inout) :: umatrix(norb,norb,norb,norb)
      complex*16, intent(in) :: rotmat(norb,norb)


      integer :: iorb1,iorb2,jorb1,jorb2,korb1,korb2,lorb1,lorb2,
     $  iorb,jorb,korb,lorb 
      complex*16 :: umatrix_temp(norb,norb,norb,norb),
     $  rotmat_temp(norb,norb)
      
      if (daggerflag .eq. 1) then
        rotmat_temp=dconjg(transpose(rotmat))
      else
        rotmat_temp=rotmat
      endif
      
      
      umatrix_temp=umatrix
      umatrix=0.0d0

      do lorb2=1, norb
        do korb2=1, norb
          do jorb2=1, norb
            do iorb2=1, norb            
              do iorb1=1, norb
                umatrix(iorb1,jorb2,korb2,lorb2)
     $            =umatrix(iorb1,jorb2,korb2,lorb2)
     $            +rotmat_temp(iorb1,iorb2)
     $            *umatrix_temp(iorb2,jorb2,korb2,lorb2)
              enddo
              
            enddo
          enddo
        enddo
      enddo

      umatrix_temp=umatrix
      umatrix=0.0d0

      do lorb2=1, norb
        do korb2=1, norb
          do jorb2=1, norb
            do iorb2=1, norb            
              do jorb1=1, norb            
                umatrix(iorb2,jorb1,korb2,lorb2)
     $            =umatrix(iorb2,jorb1,korb2,lorb2)
     $            +dconjg(rotmat_temp(jorb1,jorb2))
     $            *umatrix_temp(iorb2,jorb2,korb2,lorb2)
              enddo

            enddo
          enddo
        enddo
      enddo


      umatrix_temp=umatrix      
      umatrix=0.0d0
      do lorb2=1, norb
        do korb2=1, norb
          do jorb2=1, norb
            do iorb2=1, norb            
              do korb1=1, norb                        
                umatrix(iorb2,jorb2,korb1,lorb2)
     $            =umatrix(iorb2,jorb2,korb1,lorb2)
     $            +rotmat_temp(korb1,korb2)
     $            *umatrix_temp(iorb2,jorb2,korb2,lorb2)
              enddo

            enddo
          enddo
        enddo
      enddo


      umatrix_temp=umatrix      
      umatrix=0.0d0
      do lorb2=1, norb
        do korb2=1, norb
          do jorb2=1, norb
            do iorb2=1, norb            
              do lorb1=1, norb                        
                umatrix(iorb2,jorb2,korb2,lorb1)
     $            =umatrix(iorb2,jorb2,korb2,lorb1)
     $            +dconjg(rotmat_temp(lorb1,lorb2))
     $            *umatrix_temp(iorb2,jorb2,korb2,lorb2)
              enddo

            enddo
          enddo
        enddo
      enddo                  


c$$$  do lorb2=1, norb            
c$$$  do lorb1=1, norb                        
c$$$  do korb2=1, norb            
c$$$  do korb1=1, norb                        
c$$$  do jorb2=1, norb            
c$$$  do jorb1=1, norb            
c$$$  do iorb2=1, norb
c$$$  do iorb1=1, norb
c$$$  umatrix(iorb1,jorb1,korb1,lorb1)
c$$$  $                  =umatrix(iorb1,jorb1,korb1,lorb1)
c$$$  $                  +rotmat_temp(iorb1,iorb2)
c$$$  $                  *dconjg(rotmat_temp(jorb1,jorb2))
c$$$  $                  *rotmat_temp(korb1,korb2)
c$$$  $                  *dconjg(rotmat_temp(lorb1,lorb2))
c$$$  $                  *umatrix_temp(iorb2,jorb2,korb2,lorb2)
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  enddo
      
      end            


      

      subroutine cal_rotmat_mlms2jmj(ll, rotmat)
      
      implicit none
      integer, intent(in) :: ll
      complex*16, intent(out) :: rotmat(4*ll+2,4*ll+2)
      
      integer:: jj2,indms,indml,indmlms,indi,indmj,indjmj,
     $  iorb,jorb
      double precision :: jj,ulmu,fac1,fac2
      
      do indms=-1,1,2      
        do indml=-ll, ll
          indmlms=indml+ll+1+(indms+1)/2*(2*ll+1)
          do indi=-1,1,2
            jj2=2*ll+indi
            do indmj=-jj2, jj2, 2
              indjmj=(indmj+jj2)/2+1+(indi+1)/2*(jj2-1)
              ulmu=dble(indmj)/2.0d0/(dble(ll)+0.5d0)
              if (indi .eq. indms) then
                fac1=1.0d0
              else
                fac1=-1.0d0
              endif
              fac2=1.0d0
              if (indi .eq. -1 .and. indms .eq. 1) then
                fac2=-1.0d0
              endif
              if (indmj .eq. indms+indml*2) then
                rotmat(indjmj,indmlms)=1.0d0/dsqrt(2.0d0)
     $            *dsqrt(1+ulmu*fac1)*fac2
              endif
            enddo
          enddo
        enddo
      enddo


      end      


      subroutine cal_rotmat_cmplx2real(norb, rotmat)
      implicit none
      integer, intent(in) :: norb
      complex*16, intent(out) :: rotmat(norb,norb)

      integer :: iorb,jorb

      rotmat=0.0d0

      rotmat(norb/2+1,norb/2+1)=1.0d0       
      do iorb=1, norb/2
        rotmat(iorb+norb/2+1,iorb+norb/2+1)
     $    =1.0d0/dsqrt(2.0d0)*(-1)**iorb
        rotmat(iorb+norb/2+1,-iorb+norb/2+1)
     $    =1.0d0/dsqrt(2.0d0)
        rotmat(-iorb+norb/2+1,iorb+norb/2+1)
     $    =1.0d0/dsqrt(2.0d0)*dcmplx(0.0d0, 1.0d0)*(-1)**(iorb+1)
        rotmat(-iorb+norb/2+1,-iorb+norb/2+1)
     $    =1.0d0/dsqrt(2.0d0)*dcmplx(0.0d0, 1.0d0)
      enddo
      end      



      subroutine is_orthonormal_dcmplx(mata,lda,matb,ldb,matsize1,
     $  matsize2,matsize3,is_identity, string)

!     assumes that matsize3 < lda and ldb
!     multiplication of dconjg(transpose(a(1:matsize3,1:matsize1))) and b(1:matsize3,1:matsize2)
      use parallel_mod
      use manager_mod      
      
      implicit none
      
      integer,intent(in) :: matsize1,matsize2,matsize3,is_identity,
     $  lda,ldb
      complex*16,intent(in) :: mata(lda,matsize1),
     $  matb(ldb,matsize2)
      character(len=*), intent(in) :: string
      integer :: ii, jj
      double precision :: err, tempval
      complex*16 :: tempmat(matsize1,matsize2)

      call zgemm('c','n',matsize1,matsize2,matsize3,
     $  (1.d0,0.d0),mata,lda,
     &  matb,ldb,(0.d0,0.d0),
     $  tempmat,matsize1)

      if (is_identity .eq. 1) then
        tempval=1.0d0
      else
        tempval=0.0d0
      endif
      
      err=1.0d-6
      do ii=1, matsize1
        do jj=1, matsize2        
          if (ii .eq. jj) then
            if (cdabs(tempmat(ii,jj)-tempval) .gt. err) then
              if (maswrk) then
                write(iun,'(a, 2i5,2f12.6)')
     $            'orthonormalization fail: '//trim(string),
     $            ii, jj, tempmat(ii,jj)
                call flush(iun)
              endif
c              write(179,'(a, 2i5,2f12.6)')
c     $          'orthonormalization fail: '//trim(string),
c     $          ii, jj, tempmat(ii,jj)
c              call flush(179)              
              call ending
            endif
          else
            if (cdabs(tempmat(ii,jj)) .gt. err) then
              if (maswrk) then
c                write(iun,'(a, 2i5,2f12.6)')
c     $            'orthonormalization fail: '//trim(string),
c     $            ii, jj, tempmat(ii,jj)
c                call flush(iun)                              
              endif
c              write(179,'(a, 2i5,2f12.6)')
c     $          'orthonormalization fail: '//trim(string),
c     $          ii, jj, tempmat(ii,jj)
c              call flush(179)                            
              call ending              
            endif
          endif
        enddo
      enddo
      end            


      subroutine wan_disentangled_hamiltonian
     $  (nwan,num_bands,eigenvalues_vec,vmat,
     $  eigenvalues_dis,vmatrix_dis)

      implicit none
      
      integer, intent(in) :: nwan, num_bands
      double precision,intent(in) :: eigenvalues_vec(num_bands)
      complex*16,intent(in) :: vmat(num_bands,nwan)
      double precision, intent(out) :: eigenvalues_dis(nwan)
      complex*16, intent(out) :: vmatrix_dis(num_bands,nwan)

      integer :: iwan,jwan,ibnd

      complex*16 :: hmat(nwan,nwan)

      hmat=0.0d0
      do iwan=1, nwan
        do jwan=1, nwan
          do ibnd=1, num_bands
            hmat(iwan,jwan)=hmat(iwan,jwan)
     $        +eigenvalues_vec(ibnd)*vmat(ibnd,jwan)
     $        *dconjg(vmat(ibnd,iwan))
          enddo
        enddo
      enddo

      eigenvalues_dis=0.0d0
      call hermitianeigen_cmplxdouble(nwan,eigenvalues_dis,hmat)

      vmatrix_dis=0.0d0
      do iwan=1, nwan
        do jwan=1, nwan        
          do ibnd=1, num_bands
            vmatrix_dis(ibnd,iwan)=vmatrix_dis(ibnd,iwan)
     $        +hmat(jwan,iwan)
     $        *vmat(ibnd,jwan)
          enddo
        enddo
      enddo

      end

      SUBROUTINE  Swap_dble(a, b)
      IMPLICIT  NONE
      double precision, INTENT(INOUT) :: a, b
      double precision                :: Temp
      
      Temp = a
      a    = b
      b    = Temp
      END SUBROUTINE  Swap_dble
      
      SUBROUTINE  Swap_int(a, b)
      IMPLICIT  NONE
      integer, INTENT(INOUT) :: a, b
      integer                :: Temp
      
      Temp = a
      a    = b
      b    = Temp
      END SUBROUTINE  Swap_int
      

      
      SUBROUTINE  Sort_dble(x, Size, orderarray)
      IMPLICIT  NONE
      INTEGER, INTENT(IN) :: Size      
      double precision, INTENT(INOUT) :: x(size)
      INTEGER, INTENT(out) :: orderarray(Size)
      INTEGER :: i,orderarray_temp(Size),Location(1)
      
      do i=1, size
        orderarray_temp(i)=i
        orderarray(i)=i
      enddo
      
      DO i = 1, Size-1          ! except for the last
        Location = minloc(x(i:size))+i-1
        CALL  Swap_dble(x(i), x(Location(1)))
        CALL  Swap_int(orderarray_temp(i), orderarray_temp(Location(1)))
      END DO
      
      DO i = 1, Size-1          ! except for the last
        Location = minloc(orderarray_temp(i:size))+i-1
        CALL  Swap_int(orderarray_temp(i), orderarray_temp(Location(1)))
        CALL  Swap_int(orderarray(i), orderarray(Location(1)))
      END DO
      END



      DOUBLE PRECISION function iFactorial(j)
      IMPLICIT NONE
      INTEGER, intent(in) :: j
      INTEGER :: i
      DOUBLE PRECISION :: x
      if (j<0) print *,
     $  "iFactorial defined only for non-negative numbers!"
      x=1
      iFactorial = x
      if (j.eq.1) return
      DO i=2,j
        x = x*i
      END DO
      iFactorial = x
      return
      end function iFactorial

      DOUBLE PRECISION function dFactorial(x)
      IMPLICIT NONE
      DOUBLE PRECISION, intent(in) :: x
      DOUBLE PRECISION, PARAMETER :: spi2 = 0.8862269254527579
      DOUBLE PRECISION :: y, r
      r=1
      y=x
      DO WHILE(y.gt.1.0)
        r= r * y
        y= y -1.
      ENDDO
      IF (abs(y-0.5).LT.1e-10) r = r*spi2
      dFactorial = r
      return
      END function dFactorial

      DOUBLE PRECISION function mone(i)
      INTEGER, intent(in) :: i
      mone = 1 - 2*MOD(abs(i),2)
      return
      end function mone

      DOUBLE PRECISION function Ddelta(j1, j2, j)
      IMPLICIT NONE
      DOUBLE PRECISION, intent(in) :: j1, j2, j
!     function calls
      DOUBLE PRECISION :: dFactorial
      Ddelta = sqrt(dFactorial(j1+j2-j)*dFactorial(j1-j2+j)
     $  *dFactorial(-j1+j2+j)/dFactorial(j1+j2+j+1))
      return
      END function Ddelta



      DOUBLE PRECISION function f3j_int(j1, m1, j2, m2, j3, m3)
c$$$  wigner 3j symbol 
      IMPLICIT NONE
      integer, intent(in) :: j1, j2, j3, m1, m2, m3
      INTEGER            :: tmin, tmax, t
      DOUBLE PRECISION             :: sum, v1, v2, dn,
     $  j1d,j2d,j3d,m1d,m2d,m3d
!     function calls
      DOUBLE PRECISION             :: dFactorial
      DOUBLE PRECISION             :: iFactorial
      DOUBLE PRECISION             :: Ddelta
      DOUBLE PRECISION             :: mone
      f3j_int=0
      j1d=j1*1.0d0
      j2d=j2*1.0d0
      j3d=j3*1.0d0
      m1d=m1*1.0d0
      m2d=m2*1.0d0
      m3d=m3*1.0d0
      IF (abs(m1d+m2d+m3d) .GT. 1e-10) return
      IF (abs(j1d-j2d)-1e-14.GT.j3d.OR.j3d.GT.j1d+j2d+1e-14) return
      if (abs(m1d).GT.j1d.OR.abs(m2d).GT.j2d.OR.abs(m3d).GT.j3d)
     $  return
      tmin = INT(max(max(0.0,j2d-j3d-m1d),j1d-j3d+m2d)+1e-14)
      tmax = INT(min(min(j1d+j2d-j3d,j1d-m1d),j2d+m2d)+1e-14)
      sum=0
      DO t=tmin, tmax
        v1 = dFactorial(j3d-j2d+m1d+t)*dFactorial(j3d-j1d-m2d+t)
        v2 = dFactorial(j1d+j2d-j3d-t)*dFactorial(j1d-m1d-t)
     $    *dFactorial(j2d+m2d-t)
        sum = sum + mone(t)/(iFactorial(t)*v1*v2)
      END DO
      dn = dFactorial(j1d+m1d)*dFactorial(j1d-m1d)*dFactorial(j2d+m2d)
     $  *dFactorial(j2d-m2d)*dFactorial(j3d+m3d)*dFactorial(j3d-m3d)
      f3j_int = mone(INT(j1d-j2d-m3d))*Ddelta(j1d,j2d,j3d)*sqrt(dn)*sum
      return
      END function f3j_int      




      subroutine cal_slater_integral(mat, lval, ff)
      implicit none
      integer, intent(in) :: lval
      complex*16, intent(in) :: mat(2*lval+1,2*lval+1,2*lval+1,2*lval+1)
      complex*16, intent(out) :: ff(0:3)


      integer :: i1,i2,i3,i4,kk0,i1mi2,i3mi4,kk,mi1,mi3
      double precision, external :: f3j_int
      ff=0.0d0

      do kk=0,2*lval,2                  
        do i1=-lval, lval
          do i2=-lval, lval
            do i3=-lval, lval
              do i4=-lval, lval                  

                mi1=-i1
                mi3=-i3
                kk0=kk/2
                i1mi2=i1-i2
                i3mi4=i3-i4
                ff(kk0)=ff(kk0)
     $            +(2.0d0*kk+1.0d0)/(2.0d0*lval+1.0d0)**2
     $            /f3j_int(lval,0,kk,0,lval,0)**2
     $            *mat(i1+lval+1,i2+lval+1,i3+lval+1,i4+lval+1)
     $            *f3j_int(lval,mi1,kk,i1mi2,lval,i2)
     $            *f3j_int(lval,mi3,kk,i3mi4,lval,i4)
     $            *(-1)**(i1+i4)
              enddo
            enddo
          enddo
        enddo
      enddo
      end      

      subroutine read_wann_dat
      use wannier90_mod
      use manager_mod
      use parallel_mod
      use solid_mod

      implicit none
      include 'mpif.h'
      integer :: ik,ibnd,iwan,num_bands_temp,num_wann_temp, 
     $  ierr,ndiv_temp(3),irk,j,i,indtmp,m,mp_grid(3)
      double precision :: nelectron_dble, real_lattice_temp(3,3),
     $  recip_lattice_temp(3,3)
      complex*16, allocatable :: v_matrix_temp(:,:,:)

      if(maswrk) then

        open(unit=10, file=adjustl(trim(wanpath))//'/wannier.dat',
     $    form='unformatted')
        write(iun,*)'bni path=',
     $ adjustl(trim(wanpath))//'/wannier.dat'
        read(10) real_lattice_temp(3,3)
        read(10) recip_lattice_temp(3,3)      
        read(10) num_bands      !!
        read(10) num_wann       !!
        read(10) mp_grid
        num_kpts=mp_grid(1)*mp_grid(2)*mp_grid(3) !!
        if (.not.allocated(kpt_latt)) then
          allocate(kpt_latt(3,num_kpts),stat=ierr)
          if (ierr/=0)
     $      write(iun,*) 'Error allocating kpt_latt in read_wann_chkpt'
        endif        
        
        read(10) kpt_latt       !!
        allocate(include_bands(num_bands))
        read(10) include_bands  !!
        allocate(v_matrix(num_bands,num_wann,nqdiv))      
        read(10) v_matrix       !!
        allocate(eigenvalues(num_bands,nqdiv))
        read(10) eigenvalues    !!
        read(10) nrpts          !!
        allocate(irvec(3,nrpts))      
        read(10) irvec          !!
        allocate(ndegen(nrpts))
        read(10) ndegen         !!
        if (.not.allocated(u_matrix)) then
          allocate(u_matrix(num_wann,num_wann,num_kpts),stat=ierr)
          if (ierr/=0) write(iun,*)
     $      'Error allocating u_matrix in read_wann_chkpt'
        endif
        read(10) u_matrix       !!
        read(10) have_disentangled !!       
        write(iun,*)'have_disentangled',have_disentangled
        if (have_disentangled) then
          if (.not.allocated(lwindow)) then
            allocate(lwindow(num_bands,num_kpts),stat=ierr)
            if (ierr/=0) write(iun,*)
     $        'Error allocating lwindow in read_wann_chkpt'
          endif
          read(10) lwindow      !!

!     write(iun,*)'lwindow',lwindow
          if (.not.allocated(u_matrix_opt)) then
            allocate(u_matrix_opt(num_bands,num_wann,num_kpts),
     $        stat=ierr)
            if (ierr/=0) write(iun,*)
     $        'Error allocating u_matrix_opt in read_wann_chkpt'
          endif
          read(10) u_matrix_opt !!
!     write(iun,*)'u_matrix_opt',u_matrix_opt
        endif
        
        close(10)
      endif
      
      call mpi_bcast(num_bands,1,
     $  mpi_integer,0,MPI_COMM_WORLD,mpierr)

      call mpi_bcast(num_kpts,1,
     $  mpi_integer,0,MPI_COMM_WORLD,mpierr)      
      call mpi_bcast(num_wann,1,
     $  mpi_integer,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(have_disentangled,1,
     $  mpi_logical,0,MPI_COMM_WORLD,mpierr)
      if (.not.have_disentangled) then
        if (maswrk) then
          write(iun, *) 'have_disentangled = false not doable yet'
        endif
        call ending
      endif

      if (.not. maswrk) then
        if (have_disentangled) then        
          allocate(lwindow(num_bands,num_kpts))
          lwindow=.false.
          allocate(u_matrix_opt(num_bands,num_wann,num_kpts)) !
          u_matrix_opt=0.0d0
        endif
        allocate(u_matrix(num_wann,num_wann,num_kpts)) !
        u_matrix=0.0d0
        
        allocate(v_matrix(num_bands,num_wann,num_kpts)) !
        v_matrix=0.0d0        
        
        allocate(eigenvalues(num_bands,num_kpts)) ! 
        eigenvalues=0.0d0
        
        allocate(include_bands(num_bands))
        include_bands=0
      endif

      if (have_disentangled) then              
        call mpi_bcast(lwindow,num_bands*num_kpts,
     $    mpi_logical,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(u_matrix_opt,num_bands*num_wann*num_kpts,
     $    mpi_double_complex,0,MPI_COMM_WORLD,mpierr)
      endif

      call mpi_bcast(include_bands,num_bands,
     $  mpi_integer,0,MPI_COMM_WORLD,mpierr)      

      call mpi_bcast(u_matrix,num_wann**2*num_kpts,
     $  mpi_double_complex,0,MPI_COMM_WORLD,mpierr)      
      call mpi_bcast(v_matrix,num_wann*num_bands*num_kpts,
     $  mpi_double_complex,0,MPI_COMM_WORLD,mpierr)        

      call mpi_bcast(eigenvalues,num_bands*num_kpts,
     $  mpi_double_precision,0,MPI_COMM_WORLD,mpierr)


      call mpi_bcast(nrpts,1,
     $  mpi_integer,0,MPI_COMM_WORLD,mpierr)

      if (.not. maswrk) then
        allocate(ndegen(nrpts))
        ndegen=0
        allocate(irvec(3,nrpts))
        irvec=0
        if(allocated(kpt_latt)) deallocate(kpt_latt)
        allocate(kpt_latt(3,num_kpts))
      endif      

      call mpi_bcast(ndegen,nrpts,
     $  mpi_integer,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(irvec,3*nrpts,
     $  mpi_integer,0,MPI_COMM_WORLD,mpierr) 
      call mpi_bcast(kpt_latt,3*num_kpts,
     $  mpi_double_precision,0,MPI_COMM_WORLD,mpierr)     

      end


      subroutine za_wan_for_u(ispin,lte_bnd,gks_bnd)
      
      use parallel_mod
      use manager_mod
      use wannier90_mod
      use solid_mod
      use atom_mod
      
      implicit none
      include "mpif.h"
      
      integer, intent(in) :: ispin
      complex*16, intent(out) :: lte_bnd(nfun,num_bands,nqdiv),
     $  gks_bnd(nbasmpw,num_bands,nqdiv)
      integer :: ind_k,k,k0,iq,ibnd,jbnd
      complex*16 :: z_wan_bnd(nfun,num_bands,npnt),
     $  ev_wan_bnd(nbasmpw,num_bands,npnt)

      z_wan_bnd=(0.d0,0.d0)
      ev_wan_bnd=(0.d0,0.d0)
      

      if (me3_tau .eq. 0) then
        do ind_k=1,ndim3_k(me3_k+1)
          k=n3_mpi_k(me3_k+1)+ind_k
c$$$  do ind_k=1,ndim_k(me+1)
c$$$  k=n_mpi_k(me+1)+ind_k
          do ibnd=1, num_bands
            jbnd=include_bands(ibnd)
!     if (ubi_num .eq. 4) then
            if (ubi_num .eq. 7) then
              z_wan_bnd(:,ibnd,k)=z_qp(:,jbnd,ind_k,ispin)
              ev_wan_bnd(:,ibnd,k)=ev_qp(:,jbnd,ind_k,ispin)
!     elseif (ubi_num .gt. 4) then
!     z_wan_bnd(:,ibnd,k)=z_qp_2(:,jbnd,ind_k,ispin)
!     ev_wan_bnd(:,ibnd,k)=ev_qp_2(:,jbnd,ind_k,ispin)
            elseif (ubi_num .eq. 1) then
              z_wan_bnd(:,ibnd,k)=z_bnd(:,jbnd,ind_k,ispin)
              ev_wan_bnd(:,ibnd,k)=ev_bnd(:,jbnd,ind_k,ispin)
            endif
          enddo
        enddo
      endif
      if(nproc/=1) then
        call mpi_allreduce_dcmplx
     $    (z_wan_bnd,nfun*num_bands*npnt,mpi_sum,mpi_comm_world)
        call mpi_allreduce_dcmplx
     $    (ev_wan_bnd,nbasmpw*num_bands*npnt,mpi_sum,mpi_comm_world)
      endif        
      
      lte_bnd=0.0d0
      gks_bnd=0.0d0      

      do ind_k=1,ndim_kk(me+1)  ! k vector
        k=n_mpi_kk(me+1)+ind_k
        k0=i_kref(k)
        call sym_z_0(lte_bnd(1,1,k),k,z_wan_bnd(1,1,k0),
     &    num_bands,k_group(k),pnt(1,k))
        call sym_a_2(gks_bnd(1,1,k),k,ev_wan_bnd(1,1,k0),
     &    num_bands,k_group(k),k0)
      enddo


      if(nproc/=1) then
        call mpi_allreduce_dcmplx
     $    (lte_bnd,nfun*num_bands*nqdiv,
     $    mpi_sum,mpi_comm_world)
        call mpi_allreduce_dcmplx
     $    (gks_bnd,nbasmpw*num_bands*nqdiv,
     $    mpi_sum,mpi_comm_world)
      endif


      end
      
      
