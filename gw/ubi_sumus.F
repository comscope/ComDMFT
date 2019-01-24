      subroutine ubi_sumus(it)
      use atom_mod
	use manager_mod
	use solid_mod
	implicit none
	integer, intent(in) :: it
	integer :: i,n,isort,n1,icor
	if(it<=iter_dft) then
	  ubi='dft'
	else if(it<=iter_dft+iter_hf) then
	  ubi=' hf'
	else if(it<=iter_dft+iter_hf+iter_gw) then
	  ubi=' gw'
	else if(it<=iter_dft+iter_hf+iter_gw+iter_qp) then
	  ubi=' qp'
	endif
c ---------------------------------------
	if(iter_h_ext==0) return
      i=0
      if(ubi=='dft') then
        i=iter_h_ext/10000
      else if(ubi==' hf') then
        i=mod(iter_h_ext,10000)/1000
      else if(ubi==' qp') then
        i=mod(iter_h_ext,1000)/100
      else if(ubi==' gw') then
        i=mod(iter_h_ext,100)/10
      endif
      if(i/=0.and.nspin<nspin_0) then
        nspin=2
        n=maxmt/2
        ro(1:n)=ro(1:n)/2
        ro(n+1:n+n)=ro(1:n)
        rointr(:,1)=rointr(:,1)/2
        rointr(:,2)=rointr(:,1)
        ro_core(:,:,1)=ro_core(:,:,1)/2
        ro_core(:,:,2)=ro_core(:,:,1)
        eny(:,:,:,2)=eny(:,:,:,1)
        augm_coef(:,:,:,2)=augm_coef(:,:,:,1)
        ev_bnd(:,:,:,2)=ev_bnd(:,:,:,1)
        z_bnd(:,:,:,2)=z_bnd(:,:,:,1)
        e_bnd(:,:,2)=e_bnd(:,:,1)
        n_bnd(:,2)=n_bnd(:,1)
        do isort=1,nsort
          e_core(:,2,isort)=e_core(:,1,isort)
          r_nre_core(:,2,isort)=r_nre_core(:,1,isort)
          psi_nre(:,2,isort)=psi_nre(:,1,isort)
          cor_norm(:,2,isort)=cor_norm(:,1,isort)
        enddo
        do isort=1,nsort
	    do icor=1,ncor(isort)
            n=indcor(icor,1,isort)
            n1=indcor(icor,2,isort)
            pcor(n1:n1+nrad(isort))=pcor(n:n+nrad(isort))
            pcor_old(n1:n1+nrad(isort))=pcor_old(n:n+nrad(isort))
            qcor(n1:n1+nrad(isort))=qcor(n:n+nrad(isort))
            qcor_old(n1:n1+nrad(isort))=qcor_old(n:n+nrad(isort))
          enddo
        enddo
        g_loc_0(:,:,:,2)=g_loc_0(:,:,:,1)
        g_full_00(:,:,2)=g_full_00(:,:,1)
        if(allocated(g_full_0)) g_full_0(:,:,:,2)=g_full_0(:,:,:,1)
        if(iter_max_scf>iter_dft) then
          v_h_bb_old(:,:,:,2)=v_h_bb_old(:,:,:,1)
          v_xc_bb(:,:,:,2)=v_xc_bb(:,:,:,1)
          v_xc_bb_old(:,:,:,2)=v_xc_bb_old(:,:,:,1)
        endif
        v_intr_xc(:,2)=v_intr_xc(:,1)
        n=maxmt/2
        v_mt_xc(n+1:n+n)=v_mt_xc(1:n)
        if(iter>iter_dft+iter_hf) then
          sigx_solid_k_old(:,:,:,2)=sigx_solid_k_old(:,:,:,1)
        endif
        if(iter>iter_dft+iter_hf+1) then
          sig_c_tau_old(:,:,:,:,:,2)=sig_c_tau_old(:,:,:,:,:,1)
        endif
      endif
      end
