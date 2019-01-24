      subroutine set_control
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer :: i1,ii,isort,jsort
      itest=0
      iprn=0
c ----------------------------------------------------------------      
      n_omega=n_omega_exa+n_omega_asy
      n_nu=n_nu_exa+n_nu_asy
	iter_max_scf=iter_dft+iter_hf+iter_gw+iter_qp
c ----------------------------------------------------------------
	nqdiv=ndiv(1)*ndiv(2)*ndiv(3)
	if(maswrk) then
	  i1=len_trim(allfile)
        do ii=i1+1,72
          dosfile(ii:ii)=' '
          outfile(ii:ii)=' '
          rstfile(ii:ii)=' '
        enddo
        dosfile(1:i1+4)=allfile(1:i1)//'.dos'
        outfile(1:i1+4)=allfile(1:i1)//'.out'
        rstfile(1:i1+4)=allfile(1:i1)//'.rst'
      endif
c     Output file
      iun=2
      if(maswrk) then
         open(iun,file=outfile)
         call goto_end(iun)
      endif
c ------------------------------------------------------------------
c     limlb - number of different L,m or L,i,mj
      if(irel.ne.2) then
         limlb1=(maxb+2)**2
      else if(irel.eq.2) then
         limlb1=2*(maxb+2)**2
      endif
      limlb=nrel*(maxb+1)**2
	limlpb=(maxpb+1)**2  !! for product basis
c -----------------------------------------------------------------
      allocate(stau(3,natom))
	allocate(iob(natom))
	allocate(iob2(natom))
	allocate(iob20(nsort))
      allocate(amass(nsort))
      allocate(nhsort(nsort))
	allocate(lfun(nsort))
	allocate(lfunm(nsort))
	allocate(ncor(nsort))
      allocate(ncorm(nsort))
	allocate(lc_max(nsort))
	allocate(io_lem(natom))
	allocate(io_lem0(nsort))
	allocate(io_cor(natom))
      allocate(nlctot(nsort))
      allocate(indbas(maxntle,limlb,natom))
      allocate(indfun0(2,maxntle,nrel*maxb+1,nsort))
	call afteread
      allocate(indbasa(2,maxntle,limlb,nsort))
      allocate(lm_isz(maxel,nsort))
      allocate(lf_isz(maxel,nsort))
      allocate(ie_isz(maxel,nsort))
	allocate(indcor(ncormax,nspin_0,nsort))
	indcor=0
	allocate(ind_wf(maxlfun,nsort))
	ind_wf=0
	allocate(ntle_cor(nrel_core*maxb+1,nsort))
	allocate(indcor0(7,nrel_core*maxb+1,nsort))
	allocate(nod_cor(ncormax,nsort))
	allocate(li_cor(ncormax,nsort))
	allocate(ie_val(maxlfun,nsort))
	allocate(li_val(maxlfun,nsort))
	allocate(in_val(maxlfun,nsort))
	allocate(lm_cor(ncormmax,nsort))
	allocate(lf_cor(ncormmax,nsort))
	allocate(indcore(8,nrel_core*limlb,nsort))
	allocate(e_core(ncormax,nspin_0,nsort))
	allocate(psi_nre(ncormax,nspin_0,nsort))
	allocate(r_nre_core(ncormax,nspin_0,nsort))
	allocate(cor_norm(ncormax,nspin_0,nsort))
      allocate(iat_1(nsort))
      call flapw_1
      maxplw=400000
      nfgtl=5
      maxplwro=50000
      numvol=50
! ---- Here is for AFM irel<2 symmetrization  ---------------------      
      allocate(iafm(nsort))
      iafm=0
      if(irel<2.and.nspin==2) then
        do isort=1,nsort
          if(magn_shift(isort)<1.d-3) cycle
          do jsort=isort+1,nsort
            if(abs(z(isort)-z(jsort))>1.d-3) cycle
            if(abs(magn_shift(isort)+magn_shift(jsort))<1.d-3) then
              iafm(isort)=jsort
              iafm(jsort)=isort
            endif
          enddo
        enddo
      endif
! ----------------------------------------------------------------
      we_have_w=.false.
      omega_inf=1.d6
c ---------------------------------------------------------------      	
      v_coul_set=0
      z_fit=1
      core_mt=.true.
      w_sc_gw='scf'
      w_sc_qp='scf'
      rad_solver='DFT'
      iter_rad=1
      iter_val=1
      admix_e=0.1d0
      admix_rad=0.1d0
      end