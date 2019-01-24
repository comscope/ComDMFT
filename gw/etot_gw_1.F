      subroutine etot_gw_1(ispin,k,ind_k,gx_tau,ndim)
c -------- The VALENCE components of GW ENERGY -------------------------------------
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ndim,k,ispin,ind_k
	real*8, intent(in) :: gx_tau(ndim)
	integer :: i,j
	real*8 :: const,gcom1,conv_hh,tmp,icon,free_energy_x
	real*8, allocatable :: gb(:,:)
	allocate(gb(ndim,ndim))
	const=2.d0*wgt(k)/dfloat(nspin*nrel)
	icon=(-1)**(ispin+1)
	do i=1,ndim
	  if(ubi/='dft') gx_vh_val=gx_vh_val
     &	                        -const*gx_tau(i)*v_h_bb(i,i,ind_k,ispin)
	  ebnd_x=ebnd_x-const*gx_tau(i)*e_bnd(i,k,ispin)
	  d_free_x=d_free_x+const
     &	                 *free_energy_x(e_bnd(i,k,ispin)-chem_pot)
	enddo
	if(ubi/='dft') then
	  tmp=0.d0
        do i=1,ndim
	    tmp=tmp-const*gx_tau(i)*sigx_solid_k(i,i,ind_k,ispin)
	  enddo
	  gx_sx_val=gx_sx_val+tmp
	  e_magn=e_magn+icon*tmp
	endif
	if(ubi=='dft'.or.ubi_0=='dft') then
	  gx_vxc_val=vxch-gx_vxc_cor
	endif
	if(ubi==' qp'.or.ubi==' gw'.or.ubi=='gwg') then
c -------------- G_beta --------------------------------------	
	  if(me3_tau==0) gb=0.5d0*(gc_tau(1:ndim,1:ndim,2,1)
     &	                        -gc_tau(1:ndim,1:ndim,1,1))
        if(nproc_tau/=1) call brdcst(msgdbl,gb,8*ndim**2,0,comm_pnk)
	  if(ubi_0=='dft') then
		  gc_vxc_val=gc_vxc_val
     &	            -const*conv_hh(gb,v_xc_bb(1,1,ind_k,ispin),
     &                             ndim,ndim,nbndf)
	  endif
		gc_vh_val=gc_vh_val
     &	         -const*conv_hh(gb,v_h_bb(1,1,ind_k,ispin),
     &                          ndim,ndim,nbndf)
		tmp=-const*conv_hh(gb,sigx_solid_k(1,1,ind_k,ispin),
     &                     ndim,ndim,nbndf)
	  gc_sx_val=gc_sx_val+tmp
	  if(nproc_tau/=1) call DGOP(tmp,1,'  +',comm_pnk)
	  e_magn=e_magn+icon*tmp
	  do j=1,ndim
	    ebnd_c=ebnd_c-const*gb(j,j)*e_bnd(j,k,ispin)
	  enddo
	  call free_energy_c(ispin,ind_k)
c -------- GC_SC_VAL calculation --------------------------------------
        call int_tau_sigc_gc(gcom1,sig_c_tau(1,1,1,1,ind_k,ispin),
     &                       gc_tau,ndim,nbndf,nbndf)
	  tmp=const*gcom1*2.d0
	  gc_sc_val=gc_sc_val+tmp
	  if(nproc_tau/=1) call DGOP(tmp,1,'  +',comm_pnk)
	  e_magn=e_magn+icon*tmp
c -------- GX_SC_VAL calculation --------------------------------------
        call int_tau_sigc_gx(gcom1,sig_c_tau(1,1,1,1,ind_k,ispin),
     &                       e_bnd(1,k,ispin),ndim,nbndf)
        tmp=const*gcom1*2.d0
	  gx_sc_val=gx_sc_val+tmp
	  if(nproc_tau/=1) call DGOP(tmp,1,'  +',comm_pnk)
	  e_magn=e_magn+icon*tmp
c -------- GC_SC_VAL_W calculation -------------------------------------
        call sum_maz_sig_gc(gcom1,sig_c_omega(1,1,1,1,ind_k,ispin),
     &                      gc_omega,ndim,nbndf,nbndf)
        gc_sc_val_w=gc_sc_val_w+const*gcom1*2.d0
c -------- GX_SC_VAL_W calculation -------------------------------------
        call sum_maz_sig_gx(gcom1,sig_c_omega(1,1,1,1,ind_k,ispin),
     &                      e_bnd(1,k,ispin),ndim,nbndf)
        gx_sc_val_w=gx_sc_val_w+const*gcom1*2.d0
	endif
	deallocate(gb)
      end