      subroutine u_atomic   
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use wannier90_mod
      implicit none
      include 'mpif.h'
      integer :: iorb1,iorb2,iorb3,iorb4,ik,ispin,jspin,
     $  tau1,tau2,tau3,tau4,ind_nu,i_nu,ind_k,k,i

      complex*16, allocatable :: v_tmp(:,:),tempmat2(:,:,:,:,:,:)
C$$$      integer :: ind_nu,ind_k,iatom_c,iatom,isort,jatom_c,jatom,jsort,
C$$$     $  isp,jsp,m1,m2,m3,m4,m12,m34,ii,jj,i_nu,ni,nj,kk,
C$$$     $  aa,bb,cc,dd,ee,ff,aaindex,bbindex,mdim_i,mdim_j,mdim,
C$$$     $  i1,i2,i3,i4,l_i,l_j,mdim_loc_max,ldim_loc,li_i,li_j,li,l,
C$$$     $  ispin,jspin,pp,ind_wpbl,iwpbl,li_dim_i,li_dim_j,li_dim,
C$$$     $  iwp,jwp,ind1,ind2,tmparr(6),npi,npj,lim_i,lim_j,iq,
C$$$     $  ndim_wpbl_max,iorb,jorb,ik,iat,jat
C$$$      integer, allocatable :: lval(:)
C$$$      complex*16, allocatable ::
C$$$     $  rotmat_mlms2jmj(:,:),
C$$$     $  rotmat_cmplx2real(:,:),
C$$$     $  rotmat_temp(:,:),rotmat_temp2(:,:),
C$$$     $  uw_mat0(:,:,:,:),uw_mat1(:,:,:,:,:,:,:,:),
C$$$     $  tempmat1(:,:),tempmat2(:,:)      
      if(num_orb_proj==0) return
c     ------- We calculate interactions in 4_function basis now -----

        if (maswrk) then
          write(iun,*) 'here-1'
        endif
      
      do ind_k=1,ndim3_k(me3_k+1)
        k=n3_mpi_k(me3_k+1)+ind_k
        allocate(v_tmp(nvopt(k),nvopt(k)))
        v_tmp=0.0d0
        allocate(tempmat2(num_orb_proj,num_orb_proj,
     $    num_orb_proj,num_orb_proj,nspin,nspin))
        tempmat2=0.0d0        

        do i=1,nvopt(k)
          v_tmp(i,i)=vopt(i,ind_k)
        enddo        
        
        call boson_wan_k(v_tmp,nvopt(k),ind_k,1,1,1,
     $    v_k_den(1,1,1,1,1,ind_k),
     $    tempmat2)
        v_loc=v_loc+tempmat2
        deallocate(v_tmp)
        deallocate(tempmat2)
        
      enddo
      if (maswrk) then
          write(iun,*) 'here-2'
        endif
      call mpi_allreduce_dcmplx(v_loc,
     $  num_orb_proj**4*nspin**2,mpi_sum,comm_tau)      

      if (maswrk) then
        write(iun,*) 'here-3'
      endif
      
      call mpi_allreduce_dcmplx(w_loc,
     $  num_orb_proj**4*nspin**2*ndim3_nu,mpi_sum,comm_tau)
c$$$      call mpi_allreduce_dcmplx(epsinv_loc,
c$$$     $  num_orb_proj**4*nspin**2*ndim3_nu,mpi_sum,comm_tau)
c$$$      call mpi_allreduce_dcmplx(eps_loc,
c$$$     $  num_orb_proj**4*nspin**2*ndim3_nu,mpi_sum,comm_tau)
c$$$      call mpi_allreduce_dcmplx(chi0_loc,
c$$$     $  num_orb_proj**4*nspin**2*ndim3_nu,mpi_sum,comm_tau)
c$$$      call mpi_allreduce_dcmplx(chi_loc,
c$$$  $  num_orb_proj**4*nspin**2*ndim3_nu,mpi_sum,comm_tau)

      if (maswrk) then
        write(iun,*) 'here-4'
      endif      
      if (pol_iter_for_u .gt. 0) then
c$$$        call mpi_allreduce_dcmplx(chi0_low_loc,
c$$$     $    num_orb_proj**4*nspin**2*ndim3_nu,mpi_sum,comm_tau)      
        call mpi_allreduce_dcmplx(u_loc,
     $    num_orb_proj**4*nspin**2*ndim3_nu,mpi_sum,comm_tau)
c$$$        call mpi_allreduce_dcmplx(epsinv_high_loc,
c$$$     $    num_orb_proj**4*nspin**2*ndim3_nu,mpi_sum,comm_tau)
c$$$        call mpi_allreduce_dcmplx(eps_high_loc,
c$$$     $    num_orb_proj**4*nspin**2*ndim3_nu,mpi_sum,comm_tau)
c$$$        call mpi_allreduce_dcmplx(chi0_high_loc,
c$$$     $    num_orb_proj**4*nspin**2*ndim3_nu,mpi_sum,comm_tau)
c$$$        call mpi_allreduce_dcmplx(chi_high_loc,
c$$$     $    num_orb_proj**4*nspin**2*ndim3_nu,mpi_sum,comm_tau)
c$$$      elseif (u_mode .eq. 'edmf') then
c$$$        call mpi_allreduce_dcmplx(w_wpb_local,
c$$$     $    nwpb_max**2*natom_ciw**2*ndim3_nu,mpi_sum,comm_tau)
c$$$
c$$$      if (maswrk) then
c$$$        write(iun,*) 'here-6'
c$$$      endif        
      endif

      if (maswrk) then
        write(iun,*) 'here-7'
      endif      
c$$$  
c$$$      if (u_mode .eq. 'edmf') then
c$$$c$$$  call chi0_projection
c$$$        if (maswrk) then
c$$$          write(iun,*) 'here1'
c$$$        endif
c$$$        call u_weiss
c$$$      endif        
c$$$
c$$$      call boson_rr_dens_nu(w_k_den,'w',rydberg)
c$$$      call boson_rr_dens_nu(epsinv_k_den,'epsinv',1.0d0)
c$$$      call boson_rr_dens_nu(eps_k_den,'eps',1.0d0)
c$$$      call boson_rr_dens_nu(chi0_k_den,'chi0',1.0d0/rydberg)
c$$$      call boson_rr_dens_nu(chi_k_den,'chi',1.0d0/rydberg)
c$$$      if (pol_iter_for_u .gt. 0) then      
c$$$        call boson_rr_dens_nu(chi0_low_k_den,'chi0_low',1.0d0/rydberg)      
c$$$        call boson_rr_dens_nu(u_k_den,'u',rydberg)
c$$$        call boson_rr_dens_nu(epsinv_high_k_den,'epsinv_high',1.0d0)
c$$$        call boson_rr_dens_nu(eps_high_k_den,'eps_high',1.0d0)
c$$$        call boson_rr_dens_nu(chi0_high_k_den,'chi0_high',1.0d0/rydberg)
c$$$        call boson_rr_dens_nu(chi_high_k_den,'chi_high',1.0d0/rydberg)
c$$$        call boson_rr_dens_tau(chi_high_k_den_tau,'chi_high',
c$$$     $    1.0d0/rydberg)
c$$$      endif
c$$$      call boson_rr_dens_tau(chi_k_den_tau,'chi',
c$$$     $  1.0d0/rydberg)
c$$$      call boson_rr_dens_tau(wc_k_den_tau,'wc',
c$$$     $  rydberg)      




      if(irel<=1) then
        call v_slater_l
        call v_tensor_l
        call boson_slater_l(w_loc,'w',rydberg)
c$$$        call boson_slater_l(epsinv_loc,'epsinv',1.0d0)
c$$$        call boson_slater_l(eps_loc,'eps',1.0d0)
c$$$        call boson_slater_l(chi0_loc,'chi0',1.0d0/rydberg)
c$$$        call boson_slater_l(chi_loc,'chi',1.0d0/rydberg)
        if (pol_iter_for_u .gt. 0) then              
c$$$          call boson_slater_l(chi0_low_loc,'chi0_low',1.0d0/rydberg)        
          call boson_slater_l(u_loc,'u',rydberg)
c$$$          call boson_slater_l(epsinv_high_loc,'epsinv_high',1.0d0)
c$$$          call boson_slater_l(eps_high_loc,'eps_high',1.0d0)
c$$$          call boson_slater_l(chi0_high_loc,'chi0_high',1.0d0/rydberg)
c$$$          call boson_slater_l(chi_high_loc,'chi_high',1.0d0/rydberg)
        endif

      elseif(irel.eq.2) then
        call v_slater_rel_l
        call v_tensor_rel_l        
        call boson_slater_rel_l(w_loc,'w',rydberg)
c$$$        call boson_slater_rel_l(epsinv_loc,'epsinv',1.0d0)
c$$$        call boson_slater_rel_l(eps_loc,'eps',1.0d0)
c$$$        call boson_slater_rel_l(chi0_loc,'chi0',1.0d0/rydberg)
c$$$        call boson_slater_rel_l(chi_loc,'chi',1.0d0/rydberg)
        if (pol_iter_for_u .gt. 0) then                      
c$$$          call boson_slater_rel_l(chi0_low_loc,'chi0_low',1.0d0/rydberg)        
          call boson_slater_rel_l(u_loc,'u',rydberg)
c$$$          call boson_slater_rel_l(epsinv_high_loc,'epsinv_high',1.0d0)
c$$$          call boson_slater_rel_l(eps_high_loc,'eps_high',1.0d0)
c$$$          call boson_slater_rel_l(chi0_high_loc,'chi0_high',
c$$$     $      1.0d0/rydberg)
c$$$          call boson_slater_rel_l(chi_high_loc,'chi_high',1.0d0/rydberg)
        endif
      endif
      
c     --------------------------------------------------------------------
      end


      subroutine dcmplx_matinv(mat, invmat, dim,dimmax)
      implicit none
      integer, intent(in) :: dim, dimmax
      complex*16, intent(in) :: mat(dimmax, dimmax)
      complex*16, intent(out) :: invmat(dimmax,dimmax)
      
      integer :: ipiv(dim),info
      complex*16 :: wrk(dim**2)

      invmat=mat
      call zgetrf(dim,dim,invmat,
     $  dimmax,ipiv,info)
      if (info .ne. 0) then
c        write(179, *) 'error in zgetrf', info
      endif
      call zgetri(dim,invmat,dimmax,
     $  ipiv,wrk,dim**2,info)          
      if (info .ne. 0) then
c        write(179, *) 'error in zgetri', info
      endif
      end                
      
      
      
      subroutine boson_rr_dens_nu(k_den, tempchar,fac)
      
      use wannier90_mod
      use manager_mod      
      use parallel_mod
      use solid_mod
      use units_mod            
      

      implicit none
      include 'mpif.h'
      
      double precision,intent(in) :: fac
      complex*16,intent(in) ::
     $  k_den(num_orb_proj,nspin,num_orb_proj,nspin,
     $  k_star_max,ndim3_k(me3_k+1),ndim3_nu)
      character(len=*),intent(in) :: tempchar
      integer :: ispin,iorb,jspin,jorb,ind_nu,i_nu,ind_k,k,kst,iq,ind,ik

      complex*16 :: tmp(nqdiv,0:n_nu)

      

      open(unit=991, file=trim(adjustl(tempchar))//'_rr_den_nu.dat')                  
      open(unit=992, file=trim(adjustl(tempchar))//'_onsite_den_nu.dat')                  
      write(991,'(3(a, i5))')
     $  '#  n_nu+1= ', n_nu+1,
     $  '   nspin= ', nspin,
     $  '   nqdiv= ', nqdiv
      write(992,'(2(a, i5))')
     $  '#  n_nu+1= ', n_nu+1,
     $  '   nspin= ', nspin



      do ispin=1, nspin
        do iorb=1, num_orb_proj
          do jspin=1, nspin
            do jorb=1, num_orb_proj
              
              tmp=0.0d0                
              do ind_nu=1,ndim3_nu
                i_nu=me3_tau*ndim3_nu+ind_nu-1
                do ind_k=1, ndim3_k(me3_k+1)
                  k=n3_mpi_k(me3_k+1)+ind_k
                  do kst=1,k_star(k)
                    iq=k_list(kst,k)
                    call zone1_number(pnt(:,iq),rb0,ndiv,ind)
                    tmp(ind,i_nu)
     $                =tmp(ind,i_nu)
     $                +k_den(iorb,ispin,jorb,jspin,kst,ind_k,ind_nu)
                  enddo
                enddo
              enddo

              call mpi_allreduce_dcmplx
     $          (tmp,nqdiv*(n_nu+1),mpi_sum,mpi_comm_world)
              
              call fft3(ndiv(1),ndiv(2),ndiv(3),n_nu+1,tmp,-1)

              if (maswrk) then
                do ik=1, nqdiv
                  write(991,'(4i5,10000(2f12.6,4x))')
     $              iorb,ispin,jorb,jspin,
     $              ((tmp(ik,i_nu))*fac, i_nu=0, n_nu)
                enddo
                write(992,'(4i5,10000(2f12.6,4x))')
     $            iorb,ispin,jorb,jspin,
     $            ((tmp(1,i_nu))*fac, i_nu=0, n_nu)
              endif              
            enddo
          enddo
        enddo
      enddo
        
      if (maswrk) then
        close(991)
        close(992)        
      endif
      
      end




      subroutine boson_rr_dens_tau(k_den, tempchar,fac)
      
      use wannier90_mod
      use manager_mod      
      use parallel_mod
      use solid_mod
      use units_mod            
      

      implicit none
      include 'mpif.h'
      
      double precision,intent(in) :: fac
      complex*16,intent(in) ::
     $  k_den(num_orb_proj,nspin,num_orb_proj,nspin,
     $  k_star_max,ndim3_k(me3_k+1),ndim3_tau)
      character(len=*),intent(in) :: tempchar
      integer :: ispin,iorb,jspin,jorb,ind_tau,i_tau,
     $  ind_k,k,kst,iq,ind,ik

      complex*16 :: tmp(nqdiv,0:n_tau/2)

      

      open(unit=991, file=trim(adjustl(tempchar))//'_rr_den_tau.dat')                  
      open(unit=992, file=trim(adjustl(tempchar))//
     $  '_onsite_den_tau.dat')                  

      write(991,'(3(a, i5))')
     $  '#  n_tau/2+1= ', n_tau/2+1,
     $  '   nspin= ', nspin,
     $  '   nqdiv= ', nqdiv
      write(992,'(2(a, i5))')
     $  '#  n_tau/2+1= ', n_tau/2+1,
     $  '   nspin= ', nspin



      do ispin=1, nspin
        do iorb=1, num_orb_proj
          do jspin=1, nspin
            do jorb=1, num_orb_proj
              
              tmp=0.0d0                
              do ind_tau=1,ndim3_tau
                i_tau=me3_tau*ndim3_tau+ind_tau-1
                do ind_k=1, ndim3_k(me3_k+1)
                  k=n3_mpi_k(me3_k+1)+ind_k
                  do kst=1,k_star(k)
                    iq=k_list(kst,k)
                    call zone1_number(pnt(:,iq),rb0,ndiv,ind)
                    tmp(ind,i_tau)
     $                =tmp(ind,i_tau)
     $                +k_den(iorb,ispin,jorb,jspin,kst,ind_k,ind_tau)
                  enddo
                enddo
              enddo
!! from here
              call mpi_allreduce_dcmplx
     $          (tmp,nqdiv*(n_tau/2+1),mpi_sum,mpi_comm_world)
              
              call fft3(ndiv(1),ndiv(2),ndiv(3),n_tau/2+1,tmp,-1)

              if (maswrk) then
                do ik=1, nqdiv
                  write(991,'(4i5,10000(2f12.6,4x))')
     $              iorb,ispin,jorb,jspin,
     $              ((tmp(ik,i_tau))*fac/2.0d0, i_tau=0, n_tau/2)
                enddo
                write(992,'(4i5,10000(2f12.6,4x))')
     $            iorb,ispin,jorb,jspin,
     $            ((tmp(1,i_tau))*fac/2.0d0, i_tau=0, n_tau/2)
              endif              
            enddo
          enddo
        enddo
      enddo
        
      if (maswrk) then
        close(991)
        close(992)        
      endif
      
      end


      integer function find_ind(jj, narray,array)

      use manager_mod
      use parallel_mod
      implicit none
      integer, intent(in) :: jj,narray, array(narray)
      integer :: ii
      
      find_ind=0
      do ii=1, narray
        if (jj .eq. array(ii)) then
          find_ind=ii
        endif
      enddo
      if (find_ind .eq. 0) then
        if (maswrk) write(iun,*) 'problem in find_ind'
c        write(179,*) 'problem in find_ind'
        call ending_qft
      endif
      end



              
c$$$      subroutine u_onsite_dens(at_wpb, tempchar,is_inv)
c$$$
c$$$      use wannier90_mod
c$$$      use manager_mod      
c$$$      use parallel_mod
c$$$      use units_mod            
c$$$      
c$$$
c$$$      implicit none
c$$$      include 'mpif.h'
c$$$      
c$$$      integer, intent(in) :: is_inv
c$$$      complex*16,intent(in) :: at_wpb(nwpb,nwpb,ndim3_nu)
c$$$      character(len=*),intent(in) :: tempchar
c$$$
c$$$      integer :: aa,bb,dd,ee,i_nu,ind_nu,ispin,jspin,li_i,li_j
c$$$      double precision :: fac
c$$$      complex*16,allocatable:: tmp3(:,:),tmp4(:,:),u_mat(:,:,:,:,:)
c$$$
c$$$      allocate(tmp3(num_orb_proj*nspin,nwpb))
c$$$      tmp3=0.0d0
c$$$      allocate(tmp4(num_orb_proj*nspin,num_orb_proj*nspin))
c$$$      tmp4=0.0d0
c$$$
c$$$      allocate(u_mat(0:n_nu,num_orb_proj,nspin,
c$$$     $  num_orb_proj,nspin))
c$$$      
c$$$      u_mat=0.0d0
c$$$
c$$$      fac=rydberg
c$$$      if (is_inv .eq. 1) then
c$$$        fac=1.0d0/rydberg
c$$$      endif      
c$$$      
c$$$      do ind_nu=1,ndim3_nu
c$$$
c$$$        i_nu=me3_tau*ndim3_nu+ind_nu-1        
c$$$        
c$$$        tmp3=0.0d0
c$$$        tmp4=0.0d0
c$$$        call zgemm('n','n',num_orb_proj*nspin,nwpb,nwpb,
c$$$     $    (1.d0,0.d0),den_wanpb,num_orb_proj*nspin,
c$$$     $    at_wpb(1,1,ind_nu),
c$$$     &    nwpb,(0.d0,0.d0),tmp3,num_orb_proj*nspin)
c$$$        call zgemm('n','c',num_orb_proj*nspin,num_orb_proj*nspin,
c$$$     $    nwpb,(1.d0,0.d0),tmp3,num_orb_proj*nspin,
c$$$     $    den_wanpb(1,1,1),num_orb_proj*nspin,
c$$$     &    (0.d0,0.d0),tmp4,num_orb_proj*nspin)
c$$$        do li_i=1, num_orb_proj
c$$$          do ispin=1, nspin
c$$$            aa=(ispin-1)*num_orb_proj+li_i
c$$$            do li_j=1, num_orb_proj
c$$$              do jspin=1, nspin
c$$$                bb=(jspin-1)*num_orb_proj+li_j                    
c$$$                u_mat(i_nu,li_i,ispin,li_j,jspin)
c$$$     $            =u_mat(i_nu,li_i,ispin,li_j,jspin)
c$$$     $            +tmp4(aa,bb)
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$      deallocate(tmp3,tmp4)
c$$$      
c$$$      if(nproc_tau/=1) then
c$$$        call mpi_allreduce_dcmplx
c$$$     $    (u_mat,(num_orb_proj*nspin)**2*(n_nu+1),
c$$$     $    mpi_sum,comm_k)
c$$$      endif
c$$$      
c$$$      if (maswrk) then
c$$$        open(unit=991, file=trim(adjustl(tempchar))//'.dat')            
c$$$        write(991,'(3(a, i5))')
c$$$     $    '#  n_nu+1= ', n_nu+1,
c$$$     $    '  num_orb_proj= ', num_orb_proj,
c$$$     $    '   nspin= ', nspin
c$$$        do i_nu=0, n_nu
c$$$          write(991,'(1000f12.6)') w_nu(i_nu)*rydberg,
c$$$     $      (((( 
c$$$     $      dble(u_mat(i_nu,aa,bb,dd,ee))*fac,
c$$$     $      aa=1,num_orb_proj),
c$$$     $      bb=1,nspin),
c$$$     $      dd=1,num_orb_proj),
c$$$     $      ee=1,nspin)
c$$$        enddo
c$$$        close(991)
c$$$      endif
c$$$      deallocate(u_mat)
c$$$      end      
      


