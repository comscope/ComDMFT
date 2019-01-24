      subroutine prod_bas

      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod

      implicit none
      include "mpif.h"
      integer :: limax,ndim0,ii,jj,kk
      integer, allocatable :: nl(:,:),ind_l(:,:,:)

      limax=nrel*maxb+1
      ndim0=900

c     key=0 - Run to find the dimensions

      allocate(nl(limax,limax))
      allocate(ind_l(limax,limax,maxpb+1))
      nl=0
      ind_l=0

      call prod_b1(nl,ind_l,limax)
      if (maswrk) then
        write(993, *) 
        write(993, *) 'limax=', limax
        write(993, *) 'ind_l'
        do kk=1, maxpb+1
          write(993,*) 'kk=', kk
          do ii=1, limax
            write(993, '(100i5)') (ind_l(ii, jj, kk), jj=1, limax)
          enddo
        enddo
        write(993, *) 'nl'
        do ii=1, limax
          write(993, '(100i5)') (nl(ii, jj), jj=1, limax)
        enddo
        write(993, *) 
      endif

      call prod_b2(ndim0,0,nl,ind_l,limax)

      if(allocated(ind_pb)) deallocate(ind_pb)
      allocate(ind_pb(maxlfpb,nsort))
      if(allocated(lval_pb)) deallocate(lval_pb)
      allocate(lval_pb(maxlfpb,nsort))
      if(allocated(fun_pb)) deallocate(fun_pb)
      allocate(fun_pb(maxpbr))
      if(allocated(dfun_pb)) deallocate(dfun_pb)
      allocate(dfun_pb(max_e_pb,nsort))
      if(allocated(lm_pbmt)) deallocate(lm_pbmt)
      allocate(lm_pbmt(n_pbmtm,nsort))
      if(allocated(lf_pb)) deallocate(lf_pb)
      allocate(lf_pb(n_pbmtm,nsort))
      if(allocated(ind_prod0)) deallocate(ind_prod0)
      allocate(ind_prod0(max_e_pb,0:maxpb,nsort))
      if(allocated(ind_prod)) deallocate(ind_prod)
      allocate(ind_prod(max_e_pb,limlpb,nsort))
c     key=1 - Run to fill the arrays
      call prod_b2(ndim0,1,nl,ind_l,limax)
      if (maswrk) then
        write(993,*) 'prod_bas_test',sum(fun_pb)
      endif      
      if(allocated(fi0)) deallocate(fi0)
      allocate(fi0(maxlfun,maxlfun,maxlfpb,nsort,nspin,nspin))
      if(allocated(tmm_pb)) deallocate(tmm_pb)
      allocate(tmm_pb(max_e_pb,max_e_pb,0:maxpb,nsort))
      if(allocated(mi0)) deallocate(mi0)
      allocate(mi0(maxlfpb,nsort))
      fi0=0.d0
!      write(100,*)"fi0 before",fi0,"h",h,"nrad",nrad
      call get_fi0_gw
!      write(100,*)"fi0 after",fi0
      call mpi_barrier(mpi_comm_world,mpierr)

      call get_tmm
      if(allocated(fifi_rl)) deallocate(fifi_rl)
      allocate(fifi_rl(maxlfun,maxlfun,0:maxpb,nsort,nspin))
      call get_fifi_rl
c     ------- Calculate V_MI for k=1,npnt ------------------------------
      if(allocated(v_mi)) deallocate(v_mi)
      allocate(v_mi(n_pbmt,nplw_gw,ndim3_k(me3_k+1)))
      n_pbtot=n_pbmt+nplw_gw
      call v_matr
c     ------- Calculate S_TI_R0 for R=1,npnt ------------------------------
      call v_im_r_from_k(0,0)
      allocate(ndim_pbt(nproc_tau))
      allocate(n_mpi_pbt(nproc_tau))
      call size_shift_par(n_pbtot,nproc_tau,ndim_pbt,n_mpi_pbt)
      if(maswrk) write(iun,'(a22,i6,a1,i3,a1,a5,i6,a7,i6)')
     &  'Product basis: MT-part',n_pbmt,'(',nlb2,')',
     &  ' Int:',nplw_gw,' Total:',n_pbtot
      call v_prn
c     ------- Auxiliary array for Susceptibility calculations ----------
      if(allocated(mj_fpb)) deallocate(mj_fpb)
      allocate(mj_fpb(n_pbmtm,natom,n_symb))
!      call mj_fpb_prepare
c     --- Gradients of product basis for L=1 ---------------------------
c     do isort=1,nsort
c     do ie=1,ntle_pb(1,isort)
c     ir0=ind_prod0(ie,1,isort)  !! index of function
c     ir0=ind_pb(ir0,isort)      !! its 0-position in the full array
c     dfun_pb(ie,isort)=sqrt(0.75d0/pi)*fun_pb(ir0+1)/r(1,isort)
c     enddo
c     enddo
      deallocate(nl,ind_l)
      call timel('**** prod_bas finished *************')
      end
