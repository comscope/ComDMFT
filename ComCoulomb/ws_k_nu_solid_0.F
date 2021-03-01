      subroutine ws_k_nu_solid_0(ind_k,key_u)

      use atom_mod
      use etot_mod
c$$$  use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use wannier90_mod


      implicit none
      include "mpif.h"
      integer, intent(in) :: ind_k,key_u
      logical :: reg
      integer :: k,n_pbt,i,ind_nu,i_nu,info,j,mx2,ind_tau,ii,jj,kk,k0,
     $  ind_omega,i_omega,kst
      double precision :: pi2a,q,q2,pi8,dd1,dd2,check_hermit,dev,dev0,ak
      complex*16 :: cc,wint
      integer, allocatable :: ipiv(:)
      double precision, allocatable :: w_c_nu(:,:,:),p_low(:,:,:),
     &  p_low_tau(:,:,:)
c$$$  $  p_low_0(:,:),
c$$$  $  chimat_packed(:,:,:),chimat_high_packed(:,:,:),
c$$$  $  chimat_tau_packed(:,:,:),chimat_high_tau_packed(:,:,:)
      complex*16, allocatable :: v_tmp(:,:),v_p(:,:),w_s(:,:),
     $  epsinv(:,:),work(:,:),v_tmp1(:,:),v_tmp2(:,:),
     $  v_inv(:,:),v_p_lin(:,:),
     $  eps_singular(:,:),work1(:,:),tempmat(:,:),
     $  tempmat2(:,:,:,:,:,:)
c$$$  $  chimat_tau(:,:,:),chimat_high_tau(:,:,:),
c$$$  $  wcmat_tau(:,:,:)
      character*300 :: ch2
      double precision :: timeval
      double precision,external :: dtrace
      complex*16,external :: ztrace
      character*300 :: k0_ch, i_omega_ch,ispin_ch


c$$$      call timel('***** Ws_k_nu_solid_0 started ******')


      reg=.true.
      k=n3_mpi_k(me3_k+1)+ind_k
      pi2a=(pi+pi)/par
      pi8=8.d0*pi
      q2=pi2a**2*dot_product(pnt(:,k),pnt(:,k))
      ak=dot_product(pnt(:,k),pnt(:,k))
      q=sqrt(q2)
      n_pbt=nvopt(k)

c$$$      call cpu_time(timeval)              
c$$$      print '(a, 2i5, f12.6)', 'ws_k_nu_solid_0 start',
c$$$     $  me+1, k, timeval
c$$$      write(179, '(a, 2i5, f12.6)') 'ws_k_nu_solid_0 start',
c$$$     $  me+1, k, timeval
c$$$      call flush(179)

      allocate(work(n_pbt,n_pbt)) !
      work=0.0d0
      
      if(key_u==1) then
        allocate(p_low_tau(n_pbt,n_pbt,ndim3_tau)) !
        if (pol_iter_for_u .gt. 0) then
          do ii=1, n_pbt
            do jj=1, n_pbt
              do kk=1, ndim3_tau
                p_low_tau(ii,jj,kk)=p_wannier(ii,jj,kk,ind_k)
              enddo
            enddo
          enddo
        endif
        
      endif
      if (we_need_u) then
        allocate(tempmat2(num_orb_proj,num_orb_proj,
     $    num_orb_proj,num_orb_proj,nspin,nspin))
        tempmat2=0.0d0
      endif

c$$$  allocate(ipiv(n_pbt)) !
c$$$  ipiv=0
      allocate(work1(n_pbt,n_pbt)) !
      work1=0.0d0
      allocate(w_s(n_pbt,n_pbt)) !
      w_s=0.0d0
      allocate(tempmat(n_pbt,n_pbt)) !
      tempmat=0.0d0
      
      if(obietto=='   solid') then
        allocate(v_tmp(n_pbt,n_pbt))
        allocate(v_tmp1(n_pbt,n_pbt))
        allocate(v_tmp2(n_pbt,n_pbt))
c$$$  allocate(v_inv(n_pbt,n_pbt))
        v_tmp=0.0d0
        v_tmp1=0.0d0
        v_tmp2=0.0d0
c$$$  v_inv=0.0d0
        do i=1,n_pbt
          v_tmp(i,i)=vopt(i,ind_k)
c$$$  v_inv(i,i)=1.0d0/vopt(i,ind_k)
        enddo
c     call check_w_sym(k,v_tmp,n_pbt,2)
c     --- Contribution (P(tau=0)*V to PW-total energy ---------------------
c     dd1=0.d0
c     if(me3_tau==0) then
c     do i=1,n_pbt
c     dd1=dd1+p_w_solid_q_tau(i,i,1,ind_k)*v_tmp(i,i)
c     enddo
c     e_c_solid_pw=e_c_solid_pw-0.25d0*wgt(k)*dd1
c     endif
      endif	
      if(we_need_u .and. (pol_iter_for_u .gt. 0)) then
        allocate(p_low(n_pbt,n_pbt,ndim3_nu))
        p_low=0.0d0
        call from_tau_to_nu_ba(p_low,n_pbt,p_low_tau,n_pbt,n_pbt)
        do ind_nu=1,ndim3_nu
          call unpack_hermit(work,p_low(1,1,ind_nu),n_pbt,n_pbt,n_pbt,
     &      (0.d0,0.d0),(1.d0,0.d0))
!     for chi0_low
c$$$  call boson_wan_k(work,n_pbt,ind_k,1,1,1,
c$$$  $      chi0_low_k_den(1,1,1,1,1,ind_k,ind_nu),
c$$$  $      tempmat2)
c$$$  chi0_low_loc(:,:,:,:,:,:,ind_nu)
c$$$  $      =chi0_low_loc(:,:,:,:,:,:,ind_nu)+tempmat2

        enddo

        dev0=0.0d0
        dev=0.0d0        

        
c     --------------- REGULARIZATION of P_r-------------------------
        do ind_tau=1,ndim3_tau
c     ----- We first regularize P_LOW by removing from P_REST 
c     ------------------- positive eigen values ---
          
c$$$  if (p_sum_sign .ne. 0) then

c$$$  endif

          if(reg) then
            call boson_unpack_tau(work,n_pbt,p_low_tau(1,1,ind_tau),
     &        n_pbt,n_pbt)
c     write(179,*) 'regularization check: p_low_tau 2 vopt',
c     $        ind_k, ind_tau            
            call regularization(work,n_pbt,n_pbt,'-',dev0)
            dev=max(dev,dev0)
            call boson_pack_tau(1.d0,work,n_pbt,0.d0,
     &        p_low_tau(1,1,ind_tau),n_pbt,n_pbt)
            
            if(maswrk) write(iun,*)' Q; P_LOW: max eig.val ',k,dev
          endif
          p_low_tau(1:n_pbt,1:n_pbt,ind_tau)=
     &      p_w_solid_q_tau(1:n_pbt,1:n_pbt,ind_tau,ind_k)
     &      -p_low_tau(1:n_pbt,1:n_pbt,ind_tau)          
        enddo        
c     ------- Transform P_LOW from TAU to NU representation -------------
        
        p_low=0.0d0
        call from_tau_to_nu_ba(p_low,n_pbt,p_low_tau,n_pbt,n_pbt)
        deallocate(p_low_tau)
        if(ind_k==1) 
     &    call timel('***** P_Regularization finished ****')

c$$$  allocate(p_low_0(n_pbt,n_pbt))      
c$$$  p_low_0=0.0d0
c$$$  do ind_nu=1,ndim3_nu
c$$$  i_nu=me3_tau*ndim3_nu+ind_nu-1
c$$$  if (i_nu .eq. 0) then
c$$$  p_low_0=p_low(:,:,ind_nu)
c$$$  endif
c$$$  enddo
c$$$  call mpi_allreduce_dble(p_low_0,n_pbt**2,
c$$$  $    mpi_sum,comm_k)
      endif
      
      allocate(v_p(n_pbt,n_pbt)) !
      allocate(epsinv(n_pbt,n_pbt)) !
      if (u_linear) allocate(v_p_lin(n_pbt,n_pbt)) !
      
      
c     ------- Transform P_FULL from TAU to NU representation for given Q -

c$$$  if (we_need_w) then      
      allocate(w_c_nu(n_pbt,n_pbt,ndim3_nu)) !
c$$$  if (we_need_u) then      
c$$$  allocate(chimat_packed(n_pbt,n_pbt,ndim3_nu)) !
c$$$  chimat_packed=0.0d0
c$$$  allocate(chimat_tau_packed(n_pbt,n_pbt,ndim3_tau)) !
c$$$  chimat_tau_packed=0.0d0
c$$$  allocate(chimat_tau(n_pbt,n_pbt,ndim3_tau)) !
c$$$  chimat_tau=0.0d0
c$$$  if (pol_iter_for_u .gt. 0) then
c$$$  allocate(chimat_high_tau_packed(n_pbt,n_pbt,ndim3_tau)) !
c$$$  chimat_high_tau_packed=0.0d0        
c$$$  allocate(chimat_high_packed(n_pbt,n_pbt,ndim3_nu)) !
c$$$  chimat_high_packed=0.0d0        
c$$$  allocate(chimat_high_tau(n_pbt,n_pbt,ndim3_tau)) !
c$$$  chimat_high_tau=0.0d0
c$$$  endif
c$$$  endif
      call from_tau_to_nu_ba(w_c_nu,n_pbt,
     $  p_w_solid_q_tau(1,1,1,ind_k),nvopt_max,n_pbt)
      
c     -------- Store for output P ----------------------------------------
      call add_to_output_w(ind_k,k,n_pbt,w_c_nu,0) 
c     --------------------------------------------------------------------

c     --------------------------------------------------------------------
      dd1=0.0d0

c$$$  if(w_q_zero==1) allocate(eps_singular(2,n_pbt))
      do ind_nu=1,ndim3_nu
c$$$  if (we_need_w) then
        i_nu=me3_tau*ndim3_nu+ind_nu-1

c$$$        call cpu_time(timeval)        
c$$$        print '(a, 3i5,f12.6)', 'ws_k_nu_solid_0 k, nu',
c$$$     $    me+1, k, i_nu, timeval
c$$$        write(179, '(a, 3i5,f12.6)') 'ws_k_nu_solid_0 k, nu',
c$$$     $    me+1, k, i_nu, timeval        
c$$$        call flush(179)        


        
        work=0.0d0
        call unpack_hermit(work,w_c_nu(1,1,ind_nu),n_pbt,n_pbt,n_pbt,
     &    (0.d0,0.d0),(1.d0,0.d0))
c     --------------------------------------------------------------------
c$$$  if(w_q_zero==1) call eps_dielectric_prepare(k,ind_nu,n_pbt,work,
c$$$  &    eps_singular)
        if(q2<1.d-10) then
!     q=0 points
          if (maswrk) then
            write(993,*) 
            write(993, '(a, 4f20.6)') 'head of v_tmp at q=0 => '
     $        ,v_tmp(1,1),maxval(cdabs(v_tmp)), minval(cdabs(v_tmp))
            write(993,*) 
          endif
          
          if(w_q_zero==0) then
            if(i_nu==0.and.metal) then
              w_s=v_tmp
              call invers_h(n_pbt,w_s,n_pbt)
              w_s=w_s-work
              call invers_h(n_pbt,w_s,n_pbt)
!     if(maswrk) then
!     write(993,*)'v_tmp',v_tmp,'w_s',w_s,
!     $           'work',work,'n_pbt',n_pbt
!     endif
            else
              w_s=em1_gg0(1,i_nu)*v_tmp
            endif
          endif
c     ------- We enforce the hermiticity here ------------------------------
          dd1=max(dd1,check_hermit(w_s,n_pbt,n_pbt))
        else                    !!!if(k/=1) then
c     ------- V_P = V*P ------------------------------
c     --------- V_P = 1 - V*P ------------------------------------------
c$$$  v_p=-v_p
          if (.not. u_linear) then
            v_p=0.0d0
            epsinv=0.0d0
            do j=1,n_pbt
              do i=1,n_pbt
                v_p(i,j)=-vopt(i,ind_k)*work(i,j)
              enddo
            enddo            
            do i=1,n_pbt
              v_p(i,i)=(1.d0,0.d0)+v_p(i,i)
            enddo
          else
            v_p_lin=0.0d0
            do j=1,n_pbt
              do i=1,n_pbt
                v_p_lin(i,j)=vopt(i,ind_k)*work(i,j)
              enddo
            enddo                        
            do i=1,n_pbt
              v_p_lin(i,i)=(1.d0,0.d0)+v_p_lin(i,i)
            enddo
            v_p=0.0d0
            epsinv=0.0d0
            do j=1,n_pbt
              do i=1,n_pbt
                v_p(i,j)=-vopt(i,ind_k)*work(i,j)
              enddo
            enddo            
            do i=1,n_pbt
              v_p(i,i)=(1.d0,0.d0)+v_p(i,i)
            enddo            
          endif
c$$$  write(179,*) ind_nu,'v_p in ws_k_nu_solid_0',
c$$$  $      ztrace(n_pbt, v_p)
c     ---- We also invert V_P matrix here ------------------------------
          if(k<=npnt) then
            call pw_from_pb(ind_k,n_pbt,n_pbt,v_p,wint,2)
            e_output_int(i_nu,1,k)=dble(wint)
          endif
          e_output_int(i_nu,3,k)=dble(v_p(n_pbt,n_pbt))
          e_output_int(i_nu,4,k)=dble(v_p(1,1))            
          call dcmplx_matinv(v_p,epsinv,n_pbt,n_pbt)
          if(k<=npnt) then
            call pw_from_pb(ind_k,n_pbt,n_pbt,epsinv,wint,2)
            e_output_int(i_nu,2,k)=1.d0/dble(wint)
          endif            
c     --------- W_s = (1 - V*P)^(-1) * V   -----------------------------
          do j=1,n_pbt
            do i=1,n_pbt
              if (u_linear) then
                w_s(i,j)=v_p_lin(i,j)*vopt(j,ind_k)
              else
                w_s(i,j)=epsinv(i,j)*vopt(j,ind_k)
              endif
            enddo
          enddo
c$$$  if ((dsqrt(ak).lt.0.3d0).and.(i_nu.eq.0))then
c$$$  open(unit=991, file='static_coulomb')
c$$$  do i=1, n_pbt
c$$$  do j=1, n_pbt
c$$$  write(991,*) w_s(i,j)
c$$$  enddo
c$$$  enddo
c$$$  close(991)
c$$$  endif          
          
c$$$  write(179,*) ind_nu,'w_s in ws_k_nu_solid_0',
c$$$  $      ztrace(n_pbt, w_s)
c$$$  do ii=1, n_pbt
c$$$  write(179, *) v_p(ii, ii), w_s(ii,ii)
c$$$  enddo
c     ------- Contribution to eps^(-1) for q=(0,0,0) --------------------
          call pw_from_pb(ind_k,n_pbt,n_pbt,epsinv,cc,2)
          em1_gg0(k,i_nu)=dble(cc)
          call add_to_output_e1(ind_k,k,n_pbt,epsinv,ind_nu)
c$$$  if(w_q_zero==2) then
c$$$  c     ------- Contribution to the matrix eps^(-1) for q=(0,0,0) ---------
c$$$  call eps_invers_gamma_prep1(k,ind_nu,n_pbt,v_p)
c$$$  c     ------- Contribution to the regular part of W for q=(0,0,0) ---------
c$$$  call w_gamma_prep1(k,ind_nu,n_pbt,w_s)
c$$$  endif


        endif                   !!!if(k/=1) then
        if (we_need_u) then
          if (u_linear) then
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,
     $        (1.0d0,0.0d0),work,n_pbt,v_p_lin,n_pbt,
     $        (0.0d0,0.0d0),work1,n_pbt)
          else
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,
     $        (1.0d0,0.0d0),work,n_pbt,epsinv,n_pbt,
     $        (0.0d0,0.0d0),work1,n_pbt)
          endif

c$$$  call boson_wan_k(work1,n_pbt,ind_k,1,1,1,
c$$$  $      chi_k_den(1,1,1,1,1,ind_k,ind_nu),
c$$$  $      tempmat2)
c$$$  chi_loc(:,:,:,:,:,:,ind_nu)
c$$$  $      =chi_loc(:,:,:,:,:,:,ind_nu)+tempmat2

c$$$  call pack_hermit(work1,chimat_packed(1,1,ind_nu),
c$$$  $      n_pbt,n_pbt,n_pbt,0.d0,1.d0)

          
c     ----- Contribution to Etot_C from PW --------------------------------
! for w
          call boson_wan_k(w_s,n_pbt,ind_k,1,1,1,
     $      w_k_den(1,1,1,1,1,ind_k,ind_nu),
     $      tempmat2)
          w_loc(:,:,:,:,:,:,ind_nu)
     $      =w_loc(:,:,:,:,:,:,ind_nu)+tempmat2

c$$$  if (u_mode .eq. 'edmf') then
c$$$  call boson_wpb_k(w_s,n_pbt,ind_k,1,1,1,
c$$$  $        w_wpb_k(1,1,1,1,1,ind_k,ind_nu))
c$$$  
c$$$  do kst=1,k_star(k)                    
c$$$  w_wpb_local(:,:,:,:,ind_nu)
c$$$  $          =w_wpb_local(:,:,:,:,ind_nu)
c$$$  $          +w_wpb_k(:,:,:,:,kst,ind_k,ind_nu)          
c$$$  enddo
c$$$  endif
!     for epsilon_inv
c$$$  call boson_wan_k(epsinv,n_pbt,ind_k,1,1,1,
c$$$  $      epsinv_k_den(1,1,1,1,1,ind_k,ind_nu),
c$$$  $      tempmat2)
c$$$  epsinv_loc(:,:,:,:,:,:,ind_nu)
c$$$  $      =epsinv_loc(:,:,:,:,:,:,ind_nu)+tempmat2
!     for epsilon
c$$$  call boson_wan_k(v_p,n_pbt,ind_k,1,1,1,
c$$$  $      eps_k_den(1,1,1,1,1,ind_k,ind_nu),
c$$$  $      tempmat2)
c$$$  eps_loc(:,:,:,:,:,:,ind_nu)
c$$$  $      =eps_loc(:,:,:,:,:,:,ind_nu)+tempmat2
!     for polarizability
c$$$  call boson_wan_k(work,n_pbt,ind_k,1,1,1,
c$$$  $      chi0_k_den(1,1,1,1,1,ind_k,ind_nu),
c$$$  $      tempmat2)
c$$$  chi0_loc(:,:,:,:,:,:,ind_nu)
c$$$  $      =chi0_loc(:,:,:,:,:,:,ind_nu)+tempmat2

c$$$  call boson_wpb_k(work,n_pbt,ind_k,0,0,1,
c$$$  $      chi0_wpb_k(1,1,1,1,1,ind_k,ind_nu))
c$$$  
c$$$  do kst=1,k_star(k)          
c$$$  chi0_wpb_local(:,:,:,:,ind_nu)
c$$$  $        =chi0_wpb_local(:,:,:,:,ind_nu)
c$$$  $        +chi0_wpb_k(:,:,:,:,kst,ind_k,ind_nu)
c$$$  enddo                    
          

!     for polarizability

        endif
        w_s=w_s-v_tmp
        call etot_c_from_pw(n_pbt,k,i_nu,work,w_s)
c     -------------------------------------------------------------------
        call pack_hermit(w_s,w_c_nu(1,1,ind_nu),n_pbt,n_pbt,n_pbt,
     $    0.d0,1.d0)
c$$$  endif
c     ---------------------------------------------------------------------
        if ((we_need_u) .and. (pol_iter_for_u .gt. 0)) then
          work=0.0d0
          call unpack_hermit(work,p_low(1,1,ind_nu),n_pbt,n_pbt,n_pbt,
     &      (0.d0,0.d0),(1.d0,0.d0))

c     --------------------------------------------------------------------
          if(q2<1.d-10) then

            if (u_linear) then            
c$$$  w_s=em1_gg0_low(1,i_nu)*v_tmp
              w_s=(2.0d0-1.0d0/em1_gg0(1,i_nu))*v_tmp              
            else
              w_s=w_s+v_tmp
            endif
C$$$  if(w_q_zero==0) then
C$$$  if(i_nu==0.and.metal) then
C$$$  w_s=v_tmp
C$$$  call invers_h(n_pbt,w_s,n_pbt)
C$$$  w_s=w_s-work
C$$$  call invers_h(n_pbt,w_s,n_pbt)
C$$$  else
C$$$  w_s=em1_gg0_low(1,i_nu)*v_tmp
C$$$  endif
C$$$  endif
C$$$  c     ------- We enforce the hermiticity here ------------------------------
C$$$  dd2=max(dd2,check_hermit(w_s,n_pbt,n_pbt))
          else                  !!!if(k/=1) then
c     ------- V_P = V*P ------------------------------
            if (.not. u_linear) then
              v_p=0.0d0
              epsinv=0.0d0
              do j=1,n_pbt
                do i=1,n_pbt
                  v_p(i,j)=-vopt(i,ind_k)*work(i,j)
                enddo
              enddo            
              do i=1,n_pbt
                v_p(i,i)=(1.d0,0.d0)+v_p(i,i)
              enddo
            else
              v_p_lin=0.0d0
              do j=1,n_pbt
                do i=1,n_pbt
                  v_p_lin(i,j)=vopt(i,ind_k)*work(i,j)
                enddo
              enddo                        
              do i=1,n_pbt
                v_p_lin(i,i)=(1.d0,0.d0)+v_p_lin(i,i)
              enddo
              v_p=0.0d0
              epsinv=0.0d0
              do j=1,n_pbt
                do i=1,n_pbt
                  v_p(i,j)=-vopt(i,ind_k)*work(i,j)
                enddo
              enddo            
              do i=1,n_pbt
                v_p(i,i)=(1.d0,0.d0)+v_p(i,i)
              enddo            
            endif
c$$$  write(179,*) ind_nu,'v_p in ws_k_nu_solid_0',
c$$$  $      ztrace(n_pbt, v_p)
c     ---- We also invert V_P matrix here ------------------------------
            call dcmplx_matinv(v_p,epsinv,n_pbt,n_pbt)
c     --------- W_s = (1 - V*P)^(-1) * V   -----------------------------
            do j=1,n_pbt
              do i=1,n_pbt
                if (u_linear) then
                  w_s(i,j)=v_p_lin(i,j)*vopt(j,ind_k)
                else
                  w_s(i,j)=epsinv(i,j)*vopt(j,ind_k)
                endif
              enddo
            enddo
            call pw_from_pb(ind_k,n_pbt,n_pbt,epsinv,cc,2)
            em1_gg0_low(k,i_nu)=dble(cc)            
c     ------- Contribution to eps^(-1) for q=(0,0,0) --------------------
          endif                 !!!if(k/=1) then

          if (u_linear) then
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,
     $        (1.0d0,0.0d0),work,n_pbt,v_p_lin,n_pbt,
     $        (0.0d0,0.0d0),work1,n_pbt)
          else
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,
     $        (1.0d0,0.0d0),work,n_pbt,epsinv,n_pbt,
     $        (0.0d0,0.0d0),work1,n_pbt)
          endif

!     for polarizability
c$$$  call boson_wan_k(work1,n_pbt,ind_k,1,1,1,
c$$$  $      chi_high_k_den(1,1,1,1,1,ind_k,ind_nu),
c$$$  $      tempmat2)
c$$$  chi_high_loc(:,:,:,:,:,:,ind_nu)
c$$$  $      =chi_high_loc(:,:,:,:,:,:,ind_nu)+tempmat2
c$$$  call pack_hermit(work1,chimat_high_packed(1,1,ind_nu),
c$$$  $      n_pbt,n_pbt,n_pbt,0.d0,1.d0)
          
c     ----- Contribution to Etot_C from PW --------------------------------
          

          
          call boson_wan_k(w_s,n_pbt,ind_k,1,1,1,
     $      u_k_den(1,1,1,1,1,ind_k,ind_nu),
     $      tempmat2)          
          u_loc(:,:,:,:,:,:,ind_nu)
     $      =u_loc(:,:,:,:,:,:,ind_nu)+tempmat2
c$$$  call boson_wan_k(epsinv,n_pbt,ind_k,1,1,1,
c$$$  $      epsinv_high_k_den(1,1,1,1,1,ind_k,ind_nu),
c$$$  $      tempmat2)
c$$$  epsinv_high_loc(:,:,:,:,:,:,ind_nu)
c$$$  $      =epsinv_high_loc(:,:,:,:,:,:,ind_nu)+tempmat2
!     for epsilon
c$$$  call boson_wan_k(v_p,n_pbt,ind_k,1,1,1,
c$$$  $      eps_high_k_den(1,1,1,1,1,ind_k,ind_nu),
c$$$  $      tempmat2)
c$$$  eps_high_loc(:,:,:,:,:,:,ind_nu)
c$$$  $      =eps_high_loc(:,:,:,:,:,:,ind_nu)+tempmat2
!     for polarizability
c$$$  call boson_wan_k(work,n_pbt,ind_k,1,1,1,
c$$$  $      chi0_high_k_den(1,1,1,1,1,ind_k,ind_nu),
c$$$  $      tempmat2)
c$$$  chi0_high_loc(:,:,:,:,:,:,ind_nu)
c$$$  $      =chi0_high_loc(:,:,:,:,:,:,ind_nu)+tempmat2


c     -------------------------------------------------------------------
        endif
      enddo                     !! over ind_nu

      if(we_need_u) then
        
c$$$  chimat_tau=0.0d0
c$$$  chimat_tau_packed=0.0d0
c$$$  call from_nu_to_tau_ba(chimat_packed,n_pbt,
c$$$  $    chimat_tau_packed,n_pbt,n_pbt)
c$$$  do ind_tau=1, ndim3_tau
c$$$  call boson_unpack_tau(chimat_tau(1,1,ind_tau),n_pbt,
c$$$  &      chimat_tau_packed(1,1,ind_tau),n_pbt,n_pbt)
c$$$  call boson_wan_k(chimat_tau(1,1,ind_tau),n_pbt,ind_k,1,1,1,
c$$$  $      chi_k_den_tau(1,1,1,1,1,ind_k,ind_tau),
c$$$  $      tempmat2)
c$$$  enddo

c$$$  if (pol_iter_for_u .gt. 0) then
c$$$  chimat_high_tau=0.0d0
c$$$  chimat_high_tau_packed=0.0d0
c$$$  call from_nu_to_tau_ba(chimat_high_packed,n_pbt,
c$$$  $      chimat_high_tau_packed,n_pbt,n_pbt)
c$$$  do ind_tau=1, ndim3_tau
c$$$  call boson_unpack_tau(chimat_high_tau(1,1,ind_tau),n_pbt,
c$$$  &        chimat_high_tau_packed(1,1,ind_tau),n_pbt,n_pbt)
c$$$  call boson_wan_k(chimat_high_tau(1,1,ind_tau),
c$$$  $        n_pbt,ind_k,1,1,1,
c$$$  $        chi_high_k_den_tau(1,1,1,1,1,ind_k,ind_tau),
c$$$  $      tempmat2)
c$$$  enddo
c$$$  deallocate(chimat_high_tau_packed)        
c$$$  deallocate(chimat_high_packed)        
c$$$  deallocate(chimat_high_tau)
c$$$  endif
c$$$  deallocate(chimat_packed)
c$$$  deallocate(chimat_tau_packed)
c$$$  deallocate(chimat_tau)
        deallocate(tempmat2)        

        
      endif



      
      
c$$$  if (myrevision .eq. 1) then
c$$$  if (key_u .eq. 1) then
c$$$  call mpi_allreduce_dcmplx(p0mat_nu,n_pbt**2*(n_nu+1),
c$$$  $      mpi_sum,comm_k)
c$$$  call mpi_allreduce_dcmplx(pmat_nu,n_pbt**2*(n_nu+1),
c$$$  $      mpi_sum,comm_k)
c$$$  call mpi_allreduce_dcmplx(p0mat_low_nu,n_pbt**2*(n_nu+1),
c$$$  $      mpi_sum,comm_k)
c$$$  call mpi_allreduce_dcmplx(umat_nu,n_pbt**2*(n_nu+1),
c$$$  $      mpi_sum,comm_k)
c$$$  call mpi_allreduce_dcmplx(wmat_nu,n_pbt**2*(n_nu+1),
c$$$  $      mpi_sum,comm_k)
c$$$  
c$$$  call pade_matrix(p0mat_nu,n_pbt,k,1,p0mat_omega)
c$$$  c$$$          call pade_matrix(pmat_nu,n_pbt,k,2)
c$$$  call pade_matrix(p0mat_low_nu,n_pbt,k,3,p0mat_low_omega)
c$$$  c$$$          call pade_matrix(umat_nu,n_pbt,k,4)
c$$$  c$$$          call pade_matrix(wmat_nu,n_pbt,k,5)
c$$$  
c$$$  
c$$$  deallocate(p0mat_nu)
c$$$  deallocate(pmat_nu)
c$$$  deallocate(wmat_nu)
c$$$  deallocate(p0mat_low_nu)
c$$$  deallocate(umat_nu)
c$$$  deallocate(p0mat_omega)
c$$$  deallocate(pmat_omega)
c$$$  deallocate(wmat_omega)
c$$$  deallocate(p0mat_low_omega)
c$$$  deallocate(umat_omega)
c$$$  
c$$$  endif ! key_u
c$$$  endif ! myrevision

      deallocate(v_p)
      deallocate(epsinv)
      if (u_linear) deallocate(v_p_lin)      
c$$$  if(w_q_zero==1) deallocate(eps_singular)
      if(q2<1.d-10) then
        if(nproc_tau/=1) call mpi_allreduce_dble(dd1,1,mpi_max,comm_k)
C$$$  if(nproc_tau/=1) call mpi_allreduce_dble(dd2,1,mpi_max,comm_k)
        if(maswrk) then
          write(iun,*)' Max deviation from Hermiticity in W is ',dd1
C$$$  write(iun,*)' Max deviation from Hermiticity in U is ',dd2
        endif
      endif
c$$$  if(key_u==1) deallocate(p_low, p_low_0)
      if(ind_k==1) call timel('***** W: NU_loop finished **********')
c     -------- Store for output W ----------------------------------------
      call add_to_output_w(ind_k,k,n_pbt,w_c_nu,1)
c     ------- Transform NU ---> TAU representation ---------------------

C$$$  do ind_omega=1,ndim3_omega
C$$$  i_omega=me3_tau*ndim3_omega+ind_omega-1
C$$$  if (i_omega .lt. 10) then
C$$$  write(i_omega_ch, '(i1)') i_omega
C$$$  elseif (i_omega .lt. 100) then
C$$$  write(i_omega_ch, '(i2)') i_omega
C$$$  elseif (i_omega .lt. 1000) then
C$$$  write(i_omega_ch, '(i3)') i_omega
C$$$  elseif (i_omega .lt. 10000) then
C$$$  write(i_omega_ch, '(i4)') i_omega
C$$$  endif
C$$$  open(unit=118, file='w_'//trim(k0_ch)//'_'
C$$$  &    //trim(i_omega_ch)//'_'//'.dat')
C$$$  write(118,*) '#', n_pbt
C$$$  do i=1, n_pbt
C$$$  write(118,'(10000(f12.6,2x))')
C$$$  &      (w_c_nu(i,j,ind_omega), j=1,n_pbt)
C$$$  enddo
C$$$  close(118)
C$$$  enddo            
      call from_nu_to_tau_ba(w_c_nu,n_pbt,p_w_solid_q_tau(1,1,1,ind_k),
     &  nvopt_max,n_pbt)

c$$$  if (we_need_u) then
c$$$  allocate(tempmat2(num_orb_proj,num_orb_proj,
c$$$  $    num_orb_proj,num_orb_proj,nspin,nspin))
c$$$  
c$$$  allocate(wcmat_tau(n_pbt,n_pbt,ndim3_tau)) !
c$$$  wcmat_tau=0.0d0
c$$$  
c$$$  
c$$$  do ind_tau=1, ndim3_tau
c$$$  call boson_unpack_tau(wcmat_tau(1,1,ind_tau),n_pbt,
c$$$  &      p_w_solid_q_tau(1,1,ind_tau,ind_k),n_pbt,n_pbt)
c$$$  call boson_wan_k(wcmat_tau(1,1,ind_tau),n_pbt,ind_k,1,1,1,
c$$$  $      wc_k_den_tau(1,1,1,1,1,ind_k,ind_tau),
c$$$  $      tempmat2)
c$$$  enddo
c$$$  deallocate(wcmat_tau)
c$$$  deallocate(tempmat2)
c$$$  endif
      
      if(ind_k==1) call timel('***** W: NU_to_TAU finished ********')
      deallocate(w_c_nu)
c$$$  deallocate(ipiv)
      deallocate(work,w_s,work1,tempmat)
      if(obietto=='   solid') then
        deallocate(v_tmp,v_tmp1,v_tmp2)
c$$$  deallocate(v_inv)
      endif
      end
