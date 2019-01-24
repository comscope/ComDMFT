      subroutine mixer_gw
	use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
	use units_mod
	implicit none
	integer :: ind_k,k,ind_tau,it,n
      complex*16, allocatable :: bnd_ov_spin(:,:)
      if(iter<=iter_dft+1.or.iter>iter_max_scf) return
	if(irel<=1.and.nspin==2.and.nspin_1==1) then
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    n=n_bnd(k,1)
	    allocate(bnd_ov_spin(n,n))
	    call overl_bnd_spin(bnd_ov_spin,n,ind_k)
	    call mix_up_down(n,nbndf,v_h_bb(1,1,ind_k,1),
     &	                   v_h_bb_old(1,1,ind_k,1),v_h_bb(1,1,ind_k,2),
     &	                   v_h_bb_old(1,1,ind_k,2),bnd_ov_spin)
          if(iter>iter_dft+1) then
	      call mix_up_down(n,nbndf,sigx_solid_k(1,1,ind_k,1),
     &	                     sigx_solid_k_old(1,1,ind_k,1),
     &                       sigx_solid_k(1,1,ind_k,2),
     &                       sigx_solid_k_old(1,1,ind_k,2),
     &                       bnd_ov_spin)
          endif
          if(iter>iter_dft+iter_hf+1) then
            do ind_tau=1,ndim3_tau
              do it=1,2
	          call mix_up_down(n,nbndf,
     &	            sig_c_tau(:,:,it,ind_tau,ind_k,1),
     &	            sig_c_tau_old(:,:,it,ind_tau,ind_k,1),
     &	            sig_c_tau(:,:,it,ind_tau,ind_k,2),
     &	            sig_c_tau_old(:,:,it,ind_tau,ind_k,2),bnd_ov_spin)
              enddo
            enddo
          endif
          deallocate(bnd_ov_spin)
        enddo
      else
	  v_h_bb=adm_sig*v_h_bb+(1.d0-adm_sig)*v_h_bb_old
	  if(ubi=='dif') then
	    v_xc_bb=adm_sig*v_xc_bb+(1.d0-adm_sig)*v_xc_bb_old
	  else
          if(iter>iter_dft+1) then
	      sigx_solid_k=adm_sig*sigx_solid_k
     &	                +(1.d0-adm_sig)*sigx_solid_k_old
          endif
          if(iter>iter_dft+iter_hf+1) then
	      sig_c_tau=adm_sig*sig_c_tau+(1.d0-adm_sig)*sig_c_tau_old
          endif
        endif
      endif
      end
