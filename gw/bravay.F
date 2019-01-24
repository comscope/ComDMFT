      subroutine bravay
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: i,iplw,ia,ib,ic,k,na,nb,nc
	real*8 :: v(3),pi2,fas
c     &&&&& Here for electron gas only &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      if(maswrk) then
        write(iun,*)' The vectors of translations :'
        write(iun,'(a3,3f12.7)')' A ',(rbas(i,1),i=1,3)
        write(iun,'(a3,3f12.7)')' B ',(rbas(i,2),i=1,3)
        write(iun,'(a3,3f12.7)')' C ',(rbas(i,3),i=1,3)
      endif
c	---- Get Differences t_iatom - t_jatom  ------------
	call set_pairs
	call divide_bz
c -------------------------------------------------------------
      call getplw
c	---- Get Direct Lattice Vectors encompassing the active polygon --
      call getrbs
c ---------------------------------------------------------------------      
      call neighbors_atoms
c ----------------------------------------------------------------
      allocate(rr_pnt_ind(nqdiv,natom,natom))
      allocate(rr_pnt(3,nqdiv))
      allocate(rr_sph(3,nqdiv))
      allocate(r0_pair(nqdiv,natom,natom))
      allocate(ig_pair(nqdiv,natom,natom))
      allocate(nrr_red(3,nqdiv*natom*natom))
	call get_rr
c --------------------------------------------------------------
      allocate(indplw(-maxia:maxia,-maxib:maxib,-maxic:maxic))
	indplw=0
	do iplw=1,numplw
        ia=igbs(1,iplw)
        ib=igbs(2,iplw)
        ic=igbs(3,iplw)
        indplw(ia,ib,ic)=iplw
      enddo
c --------------------------------------------------------------
	allocate(i_kref(nqdiv))
	allocate(k_group(nqdiv))
	allocate(index_k1(nqdiv))
	allocate(ip_k(nqdiv,ngroup))
	allocate(k_star(npnt))
	allocate(ind_k_star_back(nqdiv))
	call inf_qnt(npnt,i_kref,k_group,index_k1,pnt,ip_k,k_star,
     &	           k_star_max)
      call k_sym_intact
c --------------------------------------------------------------
	call k_for_band_plot	
      call getmaxb
      call timel('***** getmaxb finished *************')
      if(maswrk) then
         write(iun,*)' Number of local basis functions ',nlcloc
         write(iun,*)' Max dimension of matrices ',
     &                       nrel*nbndtm+nlcloc
         write(iun,*)' Min dimension of matrices ',
     &                      nrel*nbndtmi+nlcloc
      endif
	allocate(indgb(nbndtm,nqdiv))
	allocate(indgb_plot(nbndtm,num_k_plot))
	call get_indgb
	call indgb_e
	allocate(nplwgw(npnt))
	allocate(nplwgw_plot(num_k_plot))
	allocate(iplf_bk(nplwbas,nqdiv))
	allocate(iplf_gk(nplwdif,nqdiv))
	allocate(k_sym_0(nqdiv,ngroup))
	allocate(g_sym_0(nqdiv,ngroup))
	allocate(gbs_sym(nplwdif,ngroup))
	call setpw_gw
	call setpw_gw1
c -------- Auxiliary array for optimization ----------------------
	if(iter_max_scf/=iter_dft) then
	  na=max(ndiv(1)-1,3)
	  nb=max(ndiv(2)-1,3)
	  nc=max(ndiv(3)-1,3)
        allocate(phase_k_rr(nqdiv,-na:na,-nb:nb,-nc:nc))	
	  pi2=pi+pi
	  do ic=-nc,nc
	    do ib=-nb,nb
	      do ia=-na,na
              v(:)=ia*rbas(:,1)+ib*rbas(:,2)+ic*rbas(:,3)
              do k=1,nqdiv
	          fas=pi2*dot_product(pnt(:,k),v)
	          phase_k_rr(k,ia,ib,ic)=dcmplx(cos(fas),sin(fas))
	        enddo
	      enddo
	    enddo
	  enddo
	endif
c --------------------------------------------------------------
      call rr_mesh_coarse
      end