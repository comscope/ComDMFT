      subroutine p_r_tau_solid(key_p)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use wannier90_mod
      use vertex_mod            
      implicit none
      include "mpif.h"      
      integer, intent(in) :: key_p
      integer :: ind_tau,ispin,key_g,ii,jj,ind_k,k,kst,iq,ind,iatom_c,
     $  n_pbt,ind_b,i_pb,jatom_c,tmparr(3),iw,jw,i_tau,i1,i2,i3,
     $  k0,i,j,ind_nu,i_nu
      double precision :: dev0,dist_r,tmppnt(3),weig,num_qpoint
      double precision, allocatable :: p_r1(:,:,:),p_r2(:,:,:)
      double precision, allocatable :: w_c_nu(:,:,:)
      complex*16 :: facdcmplx
      complex*16,allocatable :: tmp0(:,:),tmp1(:,:),
     $  tmp2(:,:),tmp3(:,:),tmp4(:,:),tmp5(:,:),tmp6(:,:),
     $  ttt1(:,:,:),ttt2(:,:,:),tmp8(:,:),work(:,:)
      character*30 :: tempch_k, tempch_tau
      character*300 :: k0_ch, i_tau_ch,ispin_ch,i_nu_ch
      double precision, external :: dtrace

      if(maswrk) 
     &  call timel('** p_r_tau_solid start *************')
c      call log_eachcore('** p_r_tau_solid start *************')



      if((ubi==' qp') .or. (ubi=='dmf')) then        
        key_g=0
      else
        key_g=1
      endif

      p_w_solid_q_tau=0.d0


      if(.not.include_gamma_w) then
        num_qpoint=nqdiv-1
        facdcmplx=dcmplx(1.0d0/dble(num_qpoint), 0.0d0)
      else
        num_qpoint=nqdiv
        facdcmplx=dcmplx(1.0d0/dble(num_qpoint), 0.0d0)        
      endif
      

      do ind_tau=1,ndim3_tau
        i_tau=me3_tau*ndim3_tau+ind_tau-1        
        do ispin=1,nspin
c     
c     ---------- MT_MT Part ------------------------------------
c     
          allocate(p_r1(n_pbmtm,n_pbmtm,nrr_max))
          p_r1=0.d0
c          call log_eachcore('** P_R_MM start *****************')
!          write(100,*)"before p_mm_r"
!          write(100,*)"p_r1",p_r1,"ind_tau",ind_tau,
!     $    "ispin",ispin,"key_g",key_g,"key_p",key_p
          call p_mm_r(p_r1,1,ind_tau,ispin,0,key_g,key_p)
          if(ind_tau==1)
     &      call timel('** P_R_MM finished *****************')
c          call log_eachcore('** P_R_MM finished *****************')
c     ---------------------------------------------------------------------
          allocate(p_r2(n_pbmt,n_pbmt,ndim3_kk(me3_k+1)))
c          write(179,*)'bkk-1'
          call sym_boson_rs(p_r1,p_r2)
c          write(179,*)'bkk-2'
          deallocate(p_r1)
          call p_mm_q_from_r(p_r2,1,ind_tau,0,key_p)
c          write(179,*)'bkk-3'
          deallocate(p_r2)
          if(ind_tau==1) 
     &      call timel('** P_MM_Q_FROM_R finished **********')
c     
c     ---------- INT_MT Part ------------------------------------
c     
          allocate(p_r2(nr_red,n_pbmt,ndim3_kk(me3_k+1)))
          p_r2=0.d0
          call p_im_r(1,ind_tau,ispin,p_r2,0,key_g,key_p)
          if(ind_tau==1) 
     &      call timel('** P_R_IM finished *****************')
c     ---------------------------------------------------------------------
          call p_im_q_from_r(p_r2,1,ind_tau,0,key_p)
          deallocate(p_r2)
          if(ind_tau==1) 
     &      call timel('** P_IM_Q_FROM_R finished **********')
        enddo                   !! over ispin
c     
c     ---------- INT_INT Part ------------------------------------
c     
        call timel('** P_ii_solid start       **********')
c        call log_eachcore('** P_ii_solid start       **********')
        call p_ii_solid(1,ind_tau,0,key_g,key_p)
c        call log_eachcore('** P_ii_solid done        **********')        
        if(ind_tau==1)
     &    call timel('** P_II_SOLID finished *************')
      enddo                     !! over ind_tau
      
      p_w_solid_q_tau=p_w_solid_q_tau*2.d0/nspin/nrel
c$$$  p_w_solid_q_tau2=p_w_solid_q_tau2*2.d0/nspin/nrel



      
c$$$  do ind_k=1,ndim3_k(me3_k+1)
c$$$  k0=n3_mpi_k(me3_k+1)+ind_k
c$$$  
c$$$  if (k0 .lt. 10) then
c$$$  write(k0_ch, '(i1)') k0
c$$$  elseif (k0 .lt. 100) then
c$$$  write(k0_ch, '(i2)') k0
c$$$  elseif (k0 .lt. 1000) then
c$$$  write(k0_ch, '(i3)') k0
c$$$  elseif (k0 .lt. 10000) then
c$$$  write(k0_ch, '(i4)') k0
c$$$  endif
c$$$  
c$$$  n_pbt=n_pbmt+nplwgw(k0)        
c$$$  allocate(w_c_nu(n_pbt,n_pbt,ndim3_nu))
c$$$  w_c_nu=0.0d0
c$$$  call from_tau_to_nu_ba
c$$$  $    (w_c_nu,n_pbt,p_w_solid_q_tau2(1,1,1,ind_k),
c$$$  &    n_pbtot,n_pbt)
c$$$  
c$$$  do ind_nu=1,ndim3_nu
c$$$  i_nu=me3_tau*ndim3_nu+ind_nu-1
c$$$  
c$$$  allocate(work(n_pbt,n_pbt))
c$$$  work=0.0d0
c$$$  call unpack_hermit(work,w_c_nu(1,1,ind_nu),n_pbt,n_pbt,n_pbt,
c$$$  $      (0.d0,0.d0),(1.d0,0.d0))
c$$$  
c$$$  C$$$          if (i_nu .lt. 10) then
c$$$  C$$$            write(i_nu_ch, '(i1)') i_nu
c$$$  C$$$          elseif (i_nu .lt. 100) then
c$$$  C$$$            write(i_nu_ch, '(i2)') i_nu
c$$$  C$$$          elseif (i_nu .lt. 1000) then
c$$$  C$$$            write(i_nu_ch, '(i3)') i_nu
c$$$  C$$$          elseif (i_nu .lt. 10000) then
c$$$  C$$$            write(i_nu_ch, '(i4)') i_nu
c$$$  C$$$          endif
c$$$  
c$$$  C$$$          open(unit=118, file='p_nu_'//trim(k0_ch)//'_'
c$$$  C$$$     &      //trim(i_nu_ch)//'_'//'.dat')
c$$$  C$$$          write(118,*) '#', n_pbt
c$$$  C$$$          do i=1, n_pbt
c$$$  C$$$            write(118,'(10000(2f12.6,2x))')
c$$$  C$$$     &        (work(i,j), j=1,n_pbt)
c$$$  C$$$          enddo
c$$$  C$$$          close(118)
c$$$  deallocate(work)
c$$$  enddo
c$$$  deallocate(w_c_nu)
c$$$  
c$$$  
c$$$  do ind_tau=1, ndim3_tau
c$$$  i_tau=me3_tau*ndim3_tau+ind_tau-1
c$$$  
c$$$  C$$$          if (i_tau .lt. 10) then
c$$$  C$$$            write(i_tau_ch, '(i1)') i_tau
c$$$  C$$$          elseif (i_tau .lt. 100) then
c$$$  C$$$            write(i_tau_ch, '(i2)') i_tau
c$$$  C$$$          elseif (i_tau .lt. 1000) then
c$$$  C$$$            write(i_tau_ch, '(i3)') i_tau
c$$$  C$$$          elseif (i_tau .lt. 10000) then
c$$$  C$$$            write(i_tau_ch, '(i4)') i_tau
c$$$  C$$$          endif
c$$$  
c$$$  C$$$          open(unit=118, file='p_tau_'//trim(k0_ch)//'_'
c$$$  C$$$     &      //trim(i_tau_ch)//'_'//'.dat')
c$$$  C$$$          write(118,*) '#', n_pbtot
c$$$  C$$$          do i=1, n_pbtot
c$$$  C$$$            write(118,'(10000(f12.6,2x))')
c$$$  C$$$     &        (p_w_solid_q_tau2(i,j,ind_tau,ind_k), j=1,n_pbtot)
c$$$  C$$$          enddo
c$$$  C$$$          close(118)
c$$$  enddo      
c$$$  enddo
      call timel('** P_R_TAU_SOLID finished **********')
      end
      




