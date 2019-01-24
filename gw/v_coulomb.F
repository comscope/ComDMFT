      subroutine v_coulomb
      use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      logical :: err
      integer :: inf,nv0,ind_q,iq,nv,nco,ind_g,ig,ind_c,ic,i,nc,n_pbt,j,
     *           gbs_number,i0,j0,ind_j,iatom,isort,ndimpb,lm,l,lget,
     &           iej,iei
      real*8 :: rout,vbz,vol,cut,v(3),v_g_q,vv,pi2a,pi8,pi2,vol0,qs,rq
      complex*16 :: cc
      integer, allocatable :: last(:),nd(:),ns(:),nd0(:),ns0(:)
      real*8, allocatable :: vec(:,:),rco(:),v0(:,:),vloc(:,:,:)
      complex*16, allocatable :: s(:,:),mt_pw(:),mt_pw_0(:,:,:),
     &                           int_pw(:),v_tmp(:,:)
	pi2=pi+pi
      pi2a=pi2/par
      pi8=8.d0*pi
c ---------------------------------------------------------------------
      vol0=pi2**3/amega/nqdiv
      qs=(3.d0*vol0/4.d0/pi)**(1.d0/3.d0)
      v_gg_00=4.d0*pi*qs/vol0
c ------ Estimate the number of G-vectors and cutoff ---------------
      nv0=numplw*16
      vbz=par**3/amega   !! in reduced units
      vol=nv0*vbz
      cut=(0.75d0*vol/pi)**(1.d0/3.d0)*0.5d0
      allocate(vec(3,nv0))
	allocate(last(0:nv0/2))
	allocate(rco(nv0/2))
	allocate(nd(nproc_tau))
	allocate(ns(nproc_tau))
	allocate(nd0(nproc_tau))
	allocate(ns0(nproc_tau))
	allocate(mt_pw(n_pbmt))
	allocate(vloc(n_pbmtm,n_pbmtm,natom))
	vloc=0.d0
	do ind_q=1,ndim3_k(me3_k+1)
	  iq=n3_mpi_k(me3_k+1)+ind_q
	  n_pbt=n_pbmt+nplwgw(iq)
	  allocate(s(nplwgw(iq),nplwgw(iq)))
	  call s_overlap(1,s,nplwgw(iq),nplwgw(iq),indpw_gw(1,iq))
c ------- List of G-vectors such that |q+G|<= cut ------------------
        call vgensp(cut,nv0,gbas,pnt(1,iq),rout,nv,vec,nco,rco,last(1),
     &              inf,err)
        last(0)=0
c ----- Number of shells included in PB_PW basis ---------------------
        nc=0
	  do i=1,nco
	    if(nplwgw(iq)==last(i)) then
	      nc=i
	      goto 2
	    endif
	  enddo    
c --------------------------------------------------------------------     
2	  call size_shift_par(nco,nproc_tau,nd0,ns0)
c -------- Prepare for MT integrals ----------------------------------
        allocate(mt_pw_0(maxlfpb,nsort,nco))
        mt_pw_0=(0.d0,0.d0)
	  do ind_c=1,nd0(me3_tau+1)
	    ic=ns0(me3_tau+1)+ind_c
          call integral_pbmt_pw_0(mt_pw_0(1,1,ic),rco(ic))
	  enddo   !! over ind_c
	  if(nproc_tau/=1) call dgop(mt_pw_0,2*maxlfpb*nsort*nco,'  +',
     & 	                           comm_pnk)
c --------------------------------------------------------------------     
	  call size_shift_par(nv,nproc_tau,nd,ns)
	  allocate(int_pw(nplwgw(iq)))
	  allocate(v_tmp(n_pbt,n_pbt))
	  v_tmp=(0.d0,0.d0)
	  do ind_g=1,nd(me3_tau+1)
	    ig=ns(me3_tau+1)+ind_g
            ic=0
	    do i=1,nco
	      if(ig<=last(i)) then
	        ic=i
	        goto 1
	      endif
	    enddo
c ------------ MT integrals --------------------------------------
1         call integral_pbmt_pw(mt_pw,mt_pw_0(1,1,ic),vec(1,ig))
c ------------ INT integrals --------------------------------------
c          if(ic>nc) call integral_pbint_pw(int_pw,vec(1,ig),iq,s)
c ----------- Original PW-form of interaction ---------------------
	    v=vec(:,ig)*pi2a
	    vv=dot_product(v,v)
            v_g_q=0.0d0
	    if(v_coul_set==0) then
	      if(vv<1.d-9) vv=1.d0/v_gg_00
	      v_g_q=pi8/vv/amega
	    else if(v_coul_set==1) then
	      if(vv<1.d-9) v_g_q=pi8*(cut_rr_eff*par)**2/amega/2
	      rq=vv*cut_rr_eff*par
	      if(vv>=1.d-9) v_g_q=pi8/vv/amega*(1.d0-cos(rq))
	    endif
	    v=vec(:,ig)-pnt(:,iq)
c ------------------------- MM --------------------------------------
          do j=1,n_pbmt
            cc=v_g_q*conjg(mt_pw(j))
            do i=1,n_pbmt
              v_tmp(i,j)=v_tmp(i,j)+mt_pw(i)*cc
            enddo
          enddo
c ------------------------- IM --------------------------------------
          if(ic<=nc) then
            i0=gbs_number(v)
            i=iplf_gk(i0,iq)+n_pbmt
            do j=1,n_pbmt
              cc=v_g_q*conjg(mt_pw(j))
              v_tmp(i,j)=v_tmp(i,j)+cc
            enddo
          else
c            do j=1,n_pbmt
c              cc=v_g_q*conjg(mt_pw(j))
c              do i=n_pbmt+1,n_pbt
c                v_tmp(i,j)=v_tmp(i,j)+int_pw(i-n_pbmt)*cc
c              enddo
c            enddo
          endif
c ------------------------- MI --------------------------------------
          if(ic<=nc) then
            j0=gbs_number(v)
            j=iplf_gk(j0,iq)+n_pbmt
            do i=1,n_pbmt
              cc=mt_pw(i)*v_g_q
              v_tmp(i,j)=v_tmp(i,j)+cc
            enddo
          else
c            do j=n_pbmt+1,n_pbt
c              cc=v_g_q*conjg(int_pw(j-n_pbmt))
c              do i=1,n_pbmt
c                v_tmp(i,j)=v_tmp(i,j)+mt_pw(i)*cc
c              enddo
c            enddo
          endif
c ------------------------- II --------------------------------------
          if(ic<=nc) then
            i0=gbs_number(v)
            i=iplf_gk(i0,iq)+n_pbmt
            v_tmp(i,i)=v_tmp(i,i)+v_g_q
          else
c            do j=n_pbmt+1,n_pbt
c              cc=v_g_q*conjg(int_pw(j-n_pbmt))
c              do i=n_pbmt+1,n_pbt
c                v_tmp(i,j)=v_tmp(i,j)+int_pw(i-n_pbmt)*cc
c              enddo
c            enddo
          endif
	  enddo   !! over ind_g
	  deallocate(mt_pw_0,int_pw,s)
c ---------- Store incomplete matrix ------------------------------
	  if(nproc_tau/=1) call dgop(v_tmp,2*n_pbt**2,'  +',comm_pnk)
	  allocate(v0(n_pbt,n_pbt))
	  call pack_hermit(v_tmp,v0,n_pbt,n_pbt,n_pbt,0.d0,1.d0)
	  do ind_j=1,ndim_pbt(me3_tau+1)
	    j=n_mpi_pbt(me3_tau+1)+ind_j
	    if(j>n_pbt) cycle
	    v_coul_pb(1:n_pbt,ind_j,ind_q)=v0(:,j)
	  enddo
c ----- Add the correction of onsite interaction ------------------
        do iatom=1,natom
	    isort=is(iatom)
	    nc=n_pbmt0(isort)
          i0=iopb(iatom)-1
          do j=1,nc
            do i=1,nc
              vloc(i,j,iatom)=vloc(i,j,iatom)+wgt(iq)*v_tmp(i0+i,i0+j)
            enddo
          enddo
        enddo	  
	  deallocate(v_tmp,v0)
      enddo   !! over ind_q
c ------ Corrected local part of interaction ------------------
	if(nproc_k/=1) call dgop(vloc,natom*n_pbmtm**2,'  +',comm_pnt)
	call sym_mm_loc_boson(vloc)
	vloc=-vloc
	do iatom=1,natom
	  isort=is(iatom)
	  ndimpb=(lmpb(isort)+1)**2
	  do lm=1,ndimpb
	    l=lget(lm)
	    do iej=1,ntle_pb(l,isort)
	      j=ind_prod(iej,lm,isort)
	      do iei=1,ntle_pb(l,isort)
	        i=ind_prod(iei,lm,isort)
	        vloc(i,j,iatom)=tmm_pb(iei,iej,l,isort)+vloc(i,j,iatom)
	      enddo
	    enddo
	  enddo
	enddo
c ------ Add local part to the full interaction matrix --------------
	do ind_q=1,ndim3_k(me3_k+1)
	  iq=n3_mpi_k(me3_k+1)+ind_q
	  n_pbt=n_pbmt+nplwgw(iq)
	  allocate(v0(n_pbt,n_pbt))
	  v0=0.d0
	  do ind_j=1,ndim_pbt(me3_tau+1)
	    j=n_mpi_pbt(me3_tau+1)+ind_j
	    if(j>n_pbt) cycle
	    v0(:,j)=v_coul_pb(1:n_pbt,ind_j,ind_q)
	  enddo
	  if(nproc_tau/=1) call dgop(v0,n_pbt**2,'  +',comm_pnk)
	  allocate(v_tmp(n_pbt,n_pbt))
	  call unpack_hermit(v_tmp,v0,n_pbt,n_pbt,n_pbt,(0.d0,0.d0),
     &	                   (1.d0,0.d0))
        do iatom=1,natom
	    isort=is(iatom)
	    nc=n_pbmt0(isort)
          i0=iopb(iatom)-1
          do j=1,nc
            do i=1,nc
              v_tmp(i0+i,i0+j)=v_tmp(i0+i,i0+j)+vloc(i,j,iatom)
            enddo
          enddo
        enddo
	  call pack_hermit(v_tmp,v0,n_pbt,n_pbt,n_pbt,0.d0,1.d0)
	  do ind_j=1,ndim_pbt(me3_tau+1)
	    j=n_mpi_pbt(me3_tau+1)+ind_j
	    if(j>n_pbt) cycle
	    v_coul_pb(1:n_pbt,ind_j,ind_q)=v0(:,j)
	  enddo
	  deallocate(v0,v_tmp)
	enddo   !! over ind_q	  	
      deallocate(vec,last,rco,nd,ns,mt_pw,nd0,ns0,vloc)
      call timel('**** V_COULOMB finished ************')
      end
