      subroutine susceptibility
      use atom_mod
      use etot_mod
      use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer :: ind_ch0(2),ind_ix0(3)
c      if(irel==2) return
c -------- LDA - based responses ----------------------------------
      if(ubi=='dft'.and.chi_cmp/='000') then
        call resp_lda_prepare(chi_cmp)
        if(nspin==1) then
          if(chi_cmp(1:1)=='1') call resp_lda('___00',1,1,1,1,0)
          if(chi_cmp(3:3)=='1') call resp_lda('___ZZ',1,1,1,n_ixc_0,0)
        else if(nspin==2) then
          if(chi_cmp(1:1)=='1'.or.chi_cmp(3:3)=='1') then
            ind_ch0(1)=ind_chi(1)
            ind_ch0(2)=ind_chi(2) 
            ind_ix0(1)=1
            ind_ix0(2)=2
            ind_ix0(3)=3   
            call resp_lda('00_ZZ',2,ind_ch0,3,ind_ix0,0)
          endif
          if(chi_cmp(2:2)=='1') then
            ind_ch0(1)=ind_chi(3)
            ind_ch0(2)=ind_chi(4) 
            ind_ix0(1)=n_ixc_0
            call resp_lda('___XY',2,ind_ch0,1,ind_ix0,0)
          endif
        endif
        deallocate(thet_me,theta_int,resp_lda_tau)
	  if(nproc_k/=1) call dgop(phi_functional_dft,1,'  +',comm_pnt)
	  if(nproc_tau/=1) call dgop(phi_functional_dft,1,'  +',comm_pnk)
      endif
c ----------------------------------------------------------
      if(ubi/='dft') then
! ------- Static (nu=0) susceptibilities along Q_LINE ------------------      
        if(vrt_x_appr(2:2)/='0') then
          if(chi_cmp(1:1)/='0'.or.chi_cmp(3:3)/='0') 
     &                                  call response_static_00_zz
        endif
! ------- Dynamic (nu>=0) susceptibilities along Q_LINE ----------------      
        if(vrt_x_appr(3:3)/='0') then
          if(chi_cmp(1:1)/='0') call response_dynamic_pw(1)
          if(chi_cmp(2:2)/='0') call response_dynamic_pw(2)
          if(chi_cmp(3:3)/='0') call response_dynamic_pw(3)
        endif
      endif
      end