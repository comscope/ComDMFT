


      subroutine cal_e_imp_rho_bnse
     $  (vmat)

      use wannier90_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use manager_mod      

      implicit none
      include 'mpif.h'

      complex*16,intent(in) ::vmat(num_bands,num_orb_proj,ndim_kk(me+1))
      
      integer :: ndiv_fi(3), nqdiv_fi,
     $  ndim_kk_fi(nproc),n_mpi_kk_fi(nproc),ndim_kk_fi_max,
     $  me1,iwan,jwan,ind_wann,wanntot,irk,iband,
     $  i1,i2,i3,j1,ik,k1,k2,k3,
     $  jj,ii,kvec(3),ibnd,
     $  ndim_nrpts(nproc),n_mpi_nrpts(nproc),ind_k,tempind(1),
     $  iorb

      double precision ::
     $  weight_band(num_bands,2), tempweightmat(num_bands)
      complex*16 :: 
     $  h_k_wan_co(nqdiv),ai,h_r_wan(nrpts)

      
      character*300 :: ch

      integer,allocatable :: wanind_cor(:),  wanind_iter(:)
      double precision,allocatable :: 
     $  occ_cor(:),e_imp_cor(:),occ_iter(:),e_imp_iter(:),
     $  occ_cor_bnd(:),e_imp_cor_bnd(:),
     $  occ_iter_bnd(:),e_imp_iter_bnd(:),
     $  h_k_wan_val(:,:)
      complex*16, allocatable ::
     $  h_k_wan_vec(:,:,:)

!      write(100,*)'kpt_latt',kpt_latt,'irvec',irvec
!      call mpi_bcast(nproc_k,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)

      ai=dcmplx(0.0d0, 1.0d0)
      me1=me+1

      allocate(occ_cor(num_orb_cor))
      allocate(e_imp_cor(num_orb_cor))
      allocate(occ_cor_bnd(num_orb_cor))
      allocate(e_imp_cor_bnd(num_orb_cor))      
      allocate(wanind_cor(num_orb_cor))
      occ_cor=0.0d0
      e_imp_cor=0.0d0
      occ_cor_bnd=0.0d0
      e_imp_cor_bnd=0.0d0      
      if (num_orb_iter .gt. 0) then
        allocate(occ_iter(num_orb_iter))
        allocate(e_imp_iter(num_orb_iter))
        allocate(occ_iter_bnd(num_orb_iter))
        allocate(e_imp_iter_bnd(num_orb_iter))        
        allocate(wanind_iter(num_orb_iter))        
        occ_iter=0.0d0
        e_imp_iter=0.0d0
        occ_iter_bnd=0.0d0
        e_imp_iter_bnd=0.0d0        
      endif

      h_r_wan=0.0d0
      
      ndiv_fi=(/20,20,20/)
      nqdiv_fi=ndiv_fi(1)*ndiv_fi(2)*ndiv_fi(3)

      call size_shift_par(nqdiv_fi,nproc,ndim_kk_fi,n_mpi_kk_fi)
      call size_shift_par(nrpts,nproc,ndim_nrpts,n_mpi_nrpts)

      allocate(h_k_wan_vec(num_orb_proj,num_orb_proj,ndim_kk_fi(me+1)))
      allocate(h_k_wan_val(num_orb_proj,ndim_kk_fi(me+1)))
      h_k_wan_vec=0.0d0
      h_k_wan_val=0.0d0
      
      ndim_kk_fi_max=maxval(ndim_kk_fi)                    

      do iwan=1, num_orb_proj
        do jwan=1, num_orb_proj
          
          h_r_wan=0.0d0
          h_k_wan_co=0.0d0
          do ind_k=1, ndim_kk(me+1)
            irk=n_mpi_kk(me+1)+ind_k          
            do iband=1, num_bands                
              h_k_wan_co(irk)=h_k_wan_co(irk)+
     $          eigenvalues(iband,irk)*vmat(iband,jwan,ind_k)
     $          *dconjg(vmat(iband,iwan,ind_k))
            enddo
          enddo
          call mpi_allreduce_dcmplx
     $      (h_k_wan_co,nqdiv,mpi_sum,mpi_comm_world)          
          
!     fourier transform to r space

          do i2=1, ndim_nrpts(me+1)
            i1=n_mpi_nrpts(me+1)+i2
            do j1=1, nqdiv
              h_r_wan(i1)=h_r_wan(i1)
     $          +h_k_wan_co(j1)*
     $          cdexp(-2.0*ai*pi*
     $          sum(irvec(:,i1)*kpt_latt(:,j1)))
     $          /dble(nqdiv)
            enddo
          enddo          
          
          call mpi_allreduce_dcmplx
     $      (h_r_wan,nrpts,mpi_sum,mpi_comm_world)
          
          do ik=1, ndim_kk_fi(me1)
            irk=n_mpi_kk_fi(me1)+ik
            call indexing(nqdiv_fi,3,ndiv_fi,0,irk,kvec)
            
            do j1=1, nrpts
              h_k_wan_vec(iwan,jwan,ik)
     $          =h_k_wan_vec(iwan,jwan,ik)
     $          +h_r_wan(j1)*
     $          cdexp(2.0*ai*pi*
     $          sum(
     $          (/dble(kvec(1)-1)/dble(ndiv_fi(1)),
     $          dble(kvec(2)-1)/dble(ndiv_fi(2)),
     $          dble(kvec(3)-1)/dble(ndiv_fi(3))/)
     $          *irvec(:,j1))
     $          )/dble(ndegen(j1))
            enddo
          enddo
        enddo
      enddo

      do ik=1, ndim_kk_fi(me1)
        call hermitianeigen_cmplxdouble(num_orb_proj,h_k_wan_val(:,ik),
     $    h_k_wan_vec(:,:,ik))

        weight_band=0.0d0
        do ibnd=1, num_orb_proj
          if (num_orb_iter .gt. 0) then
            do iorb=1, num_orb_iter
              weight_band(ibnd,2)
     $          =weight_band(ibnd,2)
     $          +cdabs(h_k_wan_vec(iorb+num_orb_cor,ibnd,ik))**2                
            enddo
          endif
          
          do iorb=1, num_orb_cor
            weight_band(ibnd,1)
     $        =weight_band(ibnd,1)
     $        +cdabs(h_k_wan_vec(iorb,ibnd,ik))**2
          enddo
        enddo
        

        tempweightmat=weight_band(:,1)
        do ii=1, num_orb_cor
          tempind=maxloc(tempweightmat)
          tempweightmat(tempind(1))=0.0d0
          wanind_cor(ii)=tempind(1)
        enddo
        
        if (num_orb_iter .gt. 0) then
          tempweightmat=weight_band(:,2)
          do ii=1, num_orb_cor
            tempweightmat(wanind_cor(ii))=0.0d0
          enddo
          do ii=1, num_orb_iter
            tempind=maxloc(tempweightmat)
            tempweightmat(tempind(1))=0.0d0
            wanind_iter(ii)=tempind(1)
          enddo
        endif        

! correlated band
        do ibnd=1, num_orb_cor
          if (h_k_wan_val(wanind_cor(ibnd),ik) .lt. 0.0d0) then
            occ_cor_bnd(ibnd)=occ_cor_bnd(ibnd)+1.0d0/dble(nqdiv_fi)
            do iwan=1, num_orb_cor
              occ_cor(iwan)=occ_cor(iwan)+1.0d0/dble(nqdiv_fi)
     $          *cdabs(h_k_wan_vec(iwan,wanind_cor(ibnd),ik))**2
            enddo
          endif            
        enddo
            
        do ibnd=1, num_orb_cor
          e_imp_cor_bnd(ibnd)=e_imp_cor_bnd(ibnd)+1.0d0/dble(nqdiv_fi)
     $      *h_k_wan_val(wanind_cor(ibnd),ik)          
          do iwan=1, num_orb_cor        
            e_imp_cor(iwan)=e_imp_cor(iwan)+1.0d0/dble(nqdiv_fi)
     $        *h_k_wan_val(wanind_cor(ibnd),ik)
     $        *cdabs(h_k_wan_vec(iwan,wanind_cor(ibnd),ik))**2
          enddo
        enddo

!     itinerant band
        if (num_orb_iter .gt. 0) then
          do ibnd=1, num_orb_iter
            if (h_k_wan_val(wanind_iter(ibnd),ik) .lt. 0.0d0) then
              occ_iter_bnd(ibnd)=occ_iter_bnd(ibnd)+1.0d0/dble(nqdiv_fi)
              do iwan=1, num_orb_iter
                occ_iter(iwan)=occ_iter(iwan)+1.0d0/dble(nqdiv_fi)
     $            *cdabs(h_k_wan_vec(iwan,wanind_iter(ibnd),ik))**2
              enddo
            endif            
          enddo
          
          do ibnd=1, num_orb_iter
            e_imp_iter_bnd(ibnd)=e_imp_iter_bnd(ibnd)
     $        +1.0d0/dble(nqdiv_fi)
     $        *h_k_wan_val(wanind_iter(ibnd),ik)
            do iwan=1, num_orb_iter        
              e_imp_iter(iwan)=e_imp_iter(iwan)+1.0d0/dble(nqdiv_fi)
     $          *h_k_wan_val(wanind_iter(ibnd),ik)
     $          *cdabs(h_k_wan_vec(iwan,wanind_iter(ibnd),ik))**2
            enddo
          enddo
        endif
        
      enddo

      call mpi_allreduce_dble
     $  (occ_cor,num_orb_cor,mpi_sum,mpi_comm_world)
      call mpi_allreduce_dble
     $  (e_imp_cor,num_orb_cor,mpi_sum,mpi_comm_world)
      call mpi_allreduce_dble
     $  (occ_cor_bnd,num_orb_cor,mpi_sum,mpi_comm_world)
      call mpi_allreduce_dble
     $  (e_imp_cor_bnd,num_orb_cor,mpi_sum,mpi_comm_world)      
      if (num_orb_iter .gt. 0) then
        call mpi_allreduce_dble
     $    (occ_iter,num_orb_iter,mpi_sum,mpi_comm_world)
        call mpi_allreduce_dble
     $    (e_imp_iter,num_orb_iter,mpi_sum,mpi_comm_world)
        call mpi_allreduce_dble
     $    (occ_iter_bnd,num_orb_iter,mpi_sum,mpi_comm_world)
        call mpi_allreduce_dble
     $    (e_imp_iter_bnd,num_orb_iter,mpi_sum,mpi_comm_world)        
      endif
      
      if (maswrk) then
        write(iun,*) "e_projected and rho_projected from fine grid"//
     $    " band-selected cor"
        do iwan=1, num_orb_cor
          write(iun,*) e_imp_cor(iwan), occ_cor(iwan)
        enddo

        if (num_orb_iter .gt. 0) then        
          write(iun,*) "e_projected and rho_projected from fine grid"//
     $      " band-selected iter"
          do iwan=1, num_orb_iter
            write(iun,*) e_imp_iter(iwan), occ_iter(iwan)
          enddo
        endif                

        write(iun,*) "e_projected and rho_projected to bands"//
     $    "from fine grid"//
     $    " band-selected cor"
        do iwan=1, num_orb_cor
          write(iun,*) e_imp_cor_bnd(iwan), occ_cor_bnd(iwan)
        enddo

        if (num_orb_iter .gt. 0) then        
          write(iun,*) "e_projected and rho_projected to bands"//
     $      " from fine grid"//
     $      " band-selected iter"
          do iwan=1, num_orb_iter
            write(iun,*) e_imp_iter_bnd(iwan), occ_iter_bnd(iwan)
          enddo
        endif        


      endif

      deallocate(occ_cor)
      deallocate(e_imp_cor)
      deallocate(occ_cor_bnd)
      deallocate(e_imp_cor_bnd)      
      deallocate(wanind_cor)
      if (num_orb_iter .gt. 0) then
        deallocate(occ_iter)
        deallocate(e_imp_iter)
        deallocate(occ_iter_bnd)
        deallocate(e_imp_iter_bnd)        
        deallocate(wanind_iter)        
      endif
      deallocate(h_k_wan_vec)
      deallocate(h_k_wan_val)      

      end                        

