      subroutine g_rs_from_ks_mm(gr,gc,it,ind_tau,ispin,key,key_p)
c$$$  eq 20 of prb 85 155129
c     -----------  G_MM from G_full in (t;b-t) form ----------------
c     key = 0 - to use QP Green's function
c     key = 1 - to use regular Green's function

      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use wannier90_mod      
      implicit none
      include "mpif.h"

      integer, intent(in) :: it,ispin,ind_tau,key,key_p
      double precision, intent(out) ::
     $  gr(maxel,maxel,ndim_nrr(me3_k+1))
      complex*16, intent(out) :: gc(maxel,maxel,ndim_nrr(me3_k+1))
      integer :: k,k0,nbnd,ind,j,ia,ib,ind_k,ig,iatom,isort,i_tau,i,
     &  i_star,j1,j0,ir,ir0,jatom,jsort,i0,i1,ind_ir0,jb
      double precision :: pi2,v(3),fc,tu,de,g_diag,green0
      complex*16 :: cc
      double precision, allocatable :: gx_tau(:)
      complex*16, allocatable :: tmp(:,:),g_tmp(:,:),tmp1(:,:),
     &  tmp2(:,:),g_mm(:,:,:),tempvec(:,:)


!        write(100,*)"chem_pot",chem_pot,"chem_pot_qp",chem_pot_qp,
!     $ "e_qp",e_qp,"e_qp_dis",e_qp_dis,"e_vec_wan",e_vec_wan,
!     $ "e_wan",e_wan,"enpj_emax",enpj_emax,"enpj_emin",enpj_emin,
!     $ "g_full",g_full,"include_bands_foru",include_bands_foru,
!     $ "ind_bands_ubn",ind_bands_ubn,"io_lem",io_lem,"irel",irel,
!     $ "is",is,"k_group",k_group,"k_list",k_list,"k_star",k_star,
!     $ "lfunm",lfunm,"maxel",maxel,"n_bnd",n_bnd,
!     $ "n_low_ubn",n_low_ubn,"n_mpi_nfun",
!     $ n_mpi_nfun,"n_mpi_nrr",n_mpi_nrr,"natom",natom,"nbndf",nbndf,
!     $ "ndim_nfun",ndim_nfun,"ndim_nrr",ndim_nrr,"ndiv",ndiv,
!     $ "nfun",nfun,"nqdiv",nqdiv,"nrr_red",nrr_red,
!     $ "num_orb_cor",num_orb_cor,"pi",pi,
!     $ "pnt",pnt,"q_qp",q_qp,"q_qp_dis",q_qp_dis,
!     $ "rb0",rb0,"tau_mesh",tau_mesh,
!     $ "tshift",tshift,"ubi",ubi,
!     $ "wan_proj_qp",wan_proj_qp,"weight_ubn",weight_ubn,
!     $ "z_bnd",z_bnd

!        write(100,*)"z_bnd",z_bnd 
!        call mpi_barrier(mpi_comm_world,mpierr)
!        call ending
      allocate(g_tmp(nfun,nfun))
      allocate(tmp2(nbndf,nfun),tmp1(nbndf,nbndf))
      allocate(tmp(nfun,nfun))
      allocate(g_mm(nqdiv,nfun,nfun))

      g_mm=(0.d0,0.d0)
      pi2=pi+pi
      i_tau=me3_tau*ndim3_tau+ind_tau-1 ! find out tau to calculate in this processor ! loop is in p_r_tau_solids.f

      if(it==2) i_tau=n_tau-i_tau ! tau pair-points
      allocate(gx_tau(nbndf))
      gx_tau=0.0d0
      
      do ind_k=1,ndim3_k(me3_k+1) ! do allocated k point calculation at each processor
c$$$  start1 : get \sum_{\lamdba,\lambda'}Z_{t,L,E}^{k,\lamdba} G_{\lambda,lambda'}^{k}(\tau) Z_{t',L',E'}^{* k,\lamdba'} in eq. 20 in IBZ
        
        k0=n3_mpi_k(me3_k+1)+ind_k ! k point index in the irreducible BZ
        nbnd=n_bnd(k0,ispin)
        g_tmp=(0.d0,0.d0)
        if(ubi=='dft'.or.ubi=='dif'.or.ubi==' hf') then
          tu=tau_mesh(i_tau)
c$$$  get g_n_tau
          call g_x_tau(ispin,k0,gx_tau,tu,nbnd,chem_pot)
c$$$  basis rotation from n(energy eigenstate) to MT basis set
          do j=1,nfun
            do ib=1,nbnd
              tmp2(ib,j)=gx_tau(ib)*conjg(z_bnd(j,ib,ind_k,ispin))
            enddo
          enddo
        else 
          if(key==0) then
            if (key_p .eq. 0) then
              tmp1=(0.d0,0.d0)
c$$$  basis rotation from n(energy eigenstate) to band basis
              do i=1,nbnd
                de=e_qp(i,k0,ispin)-chem_pot_qp
                g_diag=green0(de,tau_mesh(i_tau))
                do jb=1,nbnd
                  do ib=1,nbnd
                    tmp1(ib,jb)=tmp1(ib,jb)+q_qp(ib,i,ind_k,ispin)
     $                *g_diag*conjg(q_qp(jb,i,ind_k,ispin))
                  enddo
                enddo
              enddo
            else
              if (u_mode .eq. 'bnse') then
                tmp1=(0.d0,0.d0)
c$$$  basis rotation from n(energy eigenstate) to band basis
                do i=1, nbnd
                  if (include_bands_foru
     $              (i,ind_k,ispin,it,key_p) .eq. 1) then
                    de=e_qp(i,k0,ispin)-chem_pot_qp
                    g_diag=green0(de,tau_mesh(i_tau))
                    do jb=1,nbnd
                      do ib=1,nbnd
                        tmp1(ib,jb)=tmp1(ib,jb)
     $                    +q_qp(ib,i,ind_k,ispin)*g_diag
     &                    *conjg(q_qp(jb,i,ind_k,ispin))
                      enddo
                    enddo
                  endif
                enddo              
              elseif (u_mode .eq. 'wnds') then
                tmp1=(0.d0,0.d0)
c$$$  basis rotation from n(energy eigenstate) to band basis
                do i=1, nbnd
                  if (include_bands_foru
     $              (i,ind_k,ispin,it,key_p) .eq. 1) then
                    de=e_qp_dis(i,k0,ispin)-chem_pot_qp
                    g_diag=green0(de,tau_mesh(i_tau))
                    do jb=1,nbnd
                      do ib=1,nbnd
                        tmp1(ib,jb)=tmp1(ib,jb)
     $                    +q_qp_dis(ib,i,ind_k,ispin)*g_diag
     &                    *conjg(q_qp_dis(jb,i,ind_k,ispin))
                      enddo
                    enddo
                  endif
                enddo                  
              elseif (u_mode .eq. 'enpj') then
                tmp1=(0.d0,0.d0)
c$$$  basis rotation from n(energy eigenstate) to band basis
                do i=1,nbnd
                  de=e_qp(i,k0,ispin)-chem_pot_qp
                  if ((de .gt. enpj_emin) .and.
     $              (de .lt. enpj_emax)) then
                    g_diag=green0(de,tau_mesh(i_tau))
                    do jb=1,nbnd
                      do ib=1,nbnd
                        tmp1(ib,jb)=tmp1(ib,jb)
     $                    +q_qp(ib,i,ind_k,ispin)
     $                    *g_diag
     &                    *conjg(q_qp(jb,i,ind_k,ispin))
                      enddo
                    enddo
                  endif
                enddo
              endif
            endif
c$$$  endif              
          else if(key==1) then         
            call ferm_unpack_tau(tmp1,g_full(1,1,1,ind_tau,ind_k,ispin),
     &        nbnd,nbndf,nbndf,it)
          endif
          call zgemm('n','c',nbnd,nfun,nbnd,(1.d0,0.d0),tmp1,nbndf,
     &      z_bnd(1,1,ind_k,ispin),nfun,(0.d0,0.d0),tmp2,nbndf)
        endif
      
        call zgemm('n','n',nfun,nfun,nbnd,(1.d0,0.d0),
     &    z_bnd(1,1,ind_k,ispin),nfun,tmp2,nbndf,(0.d0,0.d0),
     &    g_tmp,nfun)           ! g_tmp: green's function at given k and tau
c$$$  end1

c$$$  start2 : evaluate \sum_{\lamdba,\lambda'}Z_{t,L,E}^{k,\lamdba} G_{\lambda,lambda'}^{k}(\tau) Z_{t',L',E'}^{* k,\lamdba'} in eq. 20 in full BZ
        do i_star=1,k_star(k0)

          k=k_list(i_star,k0)   ! k point index in the full BZ
          ig=k_group(k)
          call zone1_number(pnt(:,k),rb0,ndiv,ind) ! ind: another way to index k points in full BZ
c$$$  write(179,'(a, 5i8, 3f12.6, 2x,i8)')
c$$$  $      'zone1_number out print', ind_k, k0, i_star, k, ig,
c$$$  $      pnt(:, k), ind
          
c     ------------ Index from the left ------------------------------
          do j=1,nfun
            call sym_val_left(ig,g_tmp(:,j),tmp(:,j),1)
            do iatom=1,natom
              isort=is(iatom)
              v=tshift(:,iatom,ig)
              fc=pi2*dot_product(pnt(:,k),v)
              cc=dcmplx(cos(fc),sin(fc))
              ia=io_lem(iatom)
              ib=ia+lfunm(isort)-1
              tmp(ia:ib,j)=tmp(ia:ib,j)*cc
            enddo
          enddo
c     ------------ Index from the right  ------------------------------
          do j=1,nfun
            call sym_val_left(ig,tmp(j,:),g_mm(ind,j,:),2)
            do iatom=1,natom
              isort=is(iatom)
              v=tshift(:,iatom,ig)
              fc=-pi2*dot_product(pnt(:,k),v)
              cc=dcmplx(cos(fc),sin(fc))
              ia=io_lem(iatom)
              ib=ia+lfunm(isort)-1
              g_mm(ind,j,ia:ib)=g_mm(ind,j,ia:ib)*cc
            enddo
          enddo
        enddo
      enddo
c$$$  end2

      deallocate(tmp,g_tmp,tmp2,tmp1)
      
      if(nproc_k/=1) call mpi_allreduce_dcmplx
     $  (g_mm,nfun**2*nqdiv,mpi_sum,comm_tau)      
c     ----------------------------------------------------------
      

c$$$  start2 : evaluate eq. 20 by fourier transform
      j1=ndim_nfun(me3_k+1)
      j0=n_mpi_nfun(me3_k+1)+1
      call fft3(ndiv(1),ndiv(2),ndiv(3),nfun*j1,g_mm(1,1,j0),-1)
      do ir=1,j0-1
        g_mm(:,:,ir)=(0.d0,0.d0)
      enddo
      do ir=j0+j1,nfun
        g_mm(:,:,ir)=(0.d0,0.d0)
      enddo

      if(nproc_k/=1) call mpi_allreduce_dcmplx
     $  (g_mm,nfun**2*nqdiv,mpi_sum,comm_tau)      
      g_mm=g_mm/dfloat(nqdiv)

c     ----------------------------------------------------------
      do ind_ir0=1,ndim_nrr(me3_k+1)
        ir0=n_mpi_nrr(me3_k+1)+ind_ir0
        ir=nrr_red(1,ir0)
        iatom=nrr_red(2,ir0)
        jatom=nrr_red(3,ir0)
        isort=is(iatom)
        jsort=is(jatom)
        i0=io_lem(iatom)-1
        j0=io_lem(jatom)-1
        if(irel<=1) then
          do j1=1,lfunm(jsort)
            do i1=1,lfunm(isort)
              gr(i1,j1,ind_ir0)=g_mm(ir,i0+i1,j0+j1)
            enddo
          enddo
        else if(irel==2) then
          do j1=1,lfunm(jsort)
            do i1=1,lfunm(isort)
              gc(i1,j1,ind_ir0)=g_mm(ir,i0+i1,j0+j1)
            enddo
          enddo
        endif
      enddo
      deallocate(g_mm,gx_tau)
      end
