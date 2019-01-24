      subroutine v_spectrum

      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod

      implicit none
      include "mpif.h"
      integer :: ind_k,k,n_pbt,i,ks,i_len,n,ii
      double precision, allocatable :: vop(:,:),vpr(:,:)
      complex*16, allocatable :: v_tmp(:,:),qvop(:,:,:),s(:,:),sn(:,:)

      allocate(nvopt(npnt))
      nvopt=0
      allocate(vop(n_pbtot,ndim3_k(me3_k+1)))
      vop=0.0d0
      allocate(qvop(n_pbtot,n_pbtot,ndim3_k(me3_k+1)))
      qvop=0.0d0

      do ind_k=1,ndim3_k(me3_k+1)
        k=n3_mpi_k(me3_k+1)+ind_k
        n_pbt=n_pbmt+nplwgw(k)
        allocate(v_tmp(n_pbt,n_pbt))
c        write(179,'(a, 3i7)')
c     $    'n_pbt in v_spectrum', n_pbt, n_pbmt, nplwgw(k)
c        call flush(179)
        call timel('****  start v_coul_full, ind_k  ****')
c$$$        write(179,'(a)')
c$$$     $    'start v_coul_full'
c$$$        call flush(179)
        call v_coul_full(nplwgw(k),v_tmp,ind_k) 
c$$$        write(179,'(a)')
c$$$     $    'finish v_coul_full'
c$$$        call flush(179)
c$$$        call timel('**** finish v_coul_full, ind_k  ****')
c$$$        call timel('****  start diag_v_coul, ind_k  ****')
c$$$        write(179,'(a)')
c$$$     $    'start diag_v_coul'
c$$$        call flush(179)
        call diag_v_coul(v_tmp,n_pbt,vop(1,ind_k),qvop(1,1,ind_k),
     &    nvopt(k))
c$$$        write(179,'(a)')
c$$$     $    'finish diag_v_coul'
c$$$        call flush(179)
c$$$        call timel('**** finish diag_v_coul, ind_k  ****')
        deallocate(v_tmp)
      enddo                     !! over ind_k

      if(nproc_k/=1) call mpi_allreduce_int(nvopt,npnt,mpi_sum,comm_tau)
      if (maswrk) then
        write(993, '(a)')'nvopt'
        do ii=1, npnt
          write(993, '(2i5)') ii, nvopt(ii)
        enddo
        write(993,*)
      endif
      nvopt_max=maxval(nvopt)
      allocate(vopt(nvopt_max,ndim3_k(me3_k+1)))
      allocate(qvopt(n_pbtot,nvopt_max,ndim3_k(me3_k+1)))
      allocate(qvopt_tild(n_pbtot,nvopt_max,ndim3_k(me3_k+1)))
      allocate(qvopt_bar(n_pbtot,nvopt_max,ndim3_k(me3_k+1)))
      allocate(s(nplw_rm,nplw_rm))
      allocate(sn(nvopt_max,nvopt_max))
      vopt=0.0d0
      qvopt=0.0d0
      qvopt_tild=0.0d0
      qvopt_bar=0.0d0
      s=0.0d0
      sn=0.0d0

      do ind_k=1,ndim3_k(me3_k+1)
        k=n3_mpi_k(me3_k+1)+ind_k
        do i=1,nvopt(k)
          vopt(i,ind_k)=vop(i,ind_k)
          qvopt(:,i,ind_k)=qvop(:,i,ind_k)
        enddo
c     ---------------- Q-tilda --------------------------------------
        call s_overlap(0,s,nplw_rm,nplwgw(k),indpw_gw(1,k))
        qvopt_tild(1:n_pbmt,:,ind_k)=qvopt(1:n_pbmt,:,ind_k)
        
        call zgemm('n','n',nplwgw(k),nvopt(k),nplwgw(k),(1.d0,0.d0),s,
     &    nplw_rm,qvopt(n_pbmt+1,1,ind_k),n_pbtot,(0.d0,0.d0),
     &    qvopt_tild(n_pbmt+1,1,ind_k),n_pbtot)
c     ---------------- Q-bar --------------------------------------
        n_pbt=n_pbmt+nplwgw(k)
        call n_overlap(1,ind_k,s,nplw_rm,nplwgw(k),nvopt_max,nvopt(k),
     &    sn)
        call zgemm('n','n',n_pbt,nvopt(k),nvopt(k),(1.d0,0.d0),
     &    qvopt(1,1,ind_k),n_pbtot,sn,nvopt_max,(0.d0,0.d0),
     &    qvopt_bar(1,1,ind_k),n_pbtot)
      enddo                     !! over ind_k
      deallocate(vop,qvop,s,sn)
c     ----------------- Output ----------------------------------------
      allocate(vpr(nvopt_max,n_symb))
      vpr=0.d0
      do ind_k=1,ndim3_k(me3_k+1)
        k=n3_mpi_k(me3_k+1)+ind_k
        do ks=1,n_symb
          if(ksym_in_npnt(ks)==k) vpr(1:nvopt(k),ks)=
     &      vopt(1:nvopt(k),ind_k)
        enddo
      enddo
      if(nproc_k/=1) call mpi_allreduce_dble
     $  (vpr,nvopt_max*n_symb,mpi_sum,comm_tau)
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_V_COUL_Spectrum')
        n=min(5,n_symb)
        do ks=1,nvopt_max
          write(3,'(i5,5(1x,e12.5))')ks,(vpr(ks,i),i=1,n)
        enddo
        close(3)
      endif
      deallocate(vpr)	
      end
