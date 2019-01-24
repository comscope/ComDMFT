      subroutine output_w
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      include "mpif.h"
      integer :: k,i_len,i_nu,i,lb,idiv,i_tau
      double precision :: fc,h_r1(3),hd
      double precision, allocatable :: bos(:,:)
c     --------- Collect W_MT and W_INT  ----------------------------------
      if(nproc_k/=1) then
        call mpi_allreduce_dble
     $    (w_pr_nu,(ndiv_r1+1)*(n_nu+1),mpi_sum,comm_tau)
        call mpi_allreduce_dble
     $    (w_pr_tau,(ndiv_r1+1)*(n_tau/2+1),mpi_sum,comm_tau)
        call mpi_allreduce_dble
     $    (w_output_int,3*npnt*(n_nu+1),mpi_sum,comm_tau)
        call mpi_allreduce_dble
     $    (p_output_int,3*npnt*(n_nu+1),mpi_sum,comm_tau)
        call mpi_allreduce_dble
     $    (e_output_int,4*npnt*(n_nu+1),mpi_sum,comm_tau)
        call mpi_allreduce_dble
     $    (e1_output_int,3*npnt*(n_nu+1),mpi_sum,comm_tau)
      endif
      if(nproc_tau/=1) then
        call mpi_allreduce_dble
     $    (w_output_int,3*npnt*(n_nu+1),mpi_sum,comm_k)
        call mpi_allreduce_dble
     $    (p_output_int,3*npnt*(n_nu+1),mpi_sum,comm_k)
        call mpi_allreduce_dble
     $    (e_output_int,4*npnt*(n_nu+1),mpi_sum,comm_k)
        call mpi_allreduce_dble
     $    (e1_output_int,3*npnt*(n_nu+1),mpi_sum,comm_k)
        call mpi_allreduce_dble
     $    (w_pr_nu,(ndiv_r1+1)*(n_nu+1),mpi_sum,comm_k)
        call mpi_allreduce_dble
     $    (w_pr_tau,(ndiv_r1+1)*(n_tau/2+1),mpi_sum,comm_k)
      endif
      if(.not.maswrk) return
      fc=evolt/2.d0
      i_len=len_trim(allfile)
      if(ubi==' gw') then
        open(3,file=allfile(1:i_len)//'_W_nu_q_pw.gw')
        open(4,file=allfile(1:i_len)//'_P_nu_q_pw.gw')
        open(5,file=allfile(1:i_len)//'_E_nu_q_pw.gw')
        open(6,file=allfile(1:i_len)//'_E1_nu_q_pw.gw')
      else if(ubi=='dft') then
        open(3,file=allfile(1:i_len)//'_W_nu_q_pw.dft')
        open(4,file=allfile(1:i_len)//'_P_nu_q_pw.dft')
        open(5,file=allfile(1:i_len)//'_E_nu_q_pw.dft')
        open(6,file=allfile(1:i_len)//'_E1_nu_q_pw.dft')
      else if(ubi=='dif') then
        open(3,file=allfile(1:i_len)//'_W_nu_q_pw.dif')
        open(4,file=allfile(1:i_len)//'_P_nu_q_pw.dif')
        open(5,file=allfile(1:i_len)//'_E_nu_q_pw.dif')
        open(6,file=allfile(1:i_len)//'_E1_nu_q_pw.dif')
      else if(ubi==' hf') then
        open(3,file=allfile(1:i_len)//'_W_nu_q_pw.hf')
        open(4,file=allfile(1:i_len)//'_P_nu_q_pw.hf')
        open(5,file=allfile(1:i_len)//'_E_nu_q_pw.hf')
        open(6,file=allfile(1:i_len)//'_E1_nu_q_pw.hf')
      else if(ubi==' qp') then
        open(3,file=allfile(1:i_len)//'_W_nu_q_pw.qp')
        open(4,file=allfile(1:i_len)//'_P_nu_q_pw.qp')
        open(5,file=allfile(1:i_len)//'_E_nu_q_pw.qp')
        open(6,file=allfile(1:i_len)//'_E1_nu_q_pw.qp')
      else if(ubi=='gwg') then
        open(3,file=allfile(1:i_len)//'_W_nu_q_pw.gwg')
        open(4,file=allfile(1:i_len)//'_P_nu_q_pw.gwg')
        open(5,file=allfile(1:i_len)//'_E_nu_q_pw.gwg')
        open(6,file=allfile(1:i_len)//'_E1_nu_q_pw.gwg')
      else if(ubi=='dmf') then
        open(3,file=allfile(1:i_len)//'_W_nu_q_pw.dmf')
        open(4,file=allfile(1:i_len)//'_P_nu_q_pw.dmf')
        open(5,file=allfile(1:i_len)//'_E_nu_q_pw.dmf')
        open(6,file=allfile(1:i_len)//'_E1_nu_q_pw.dmf')
      endif
      do k=1,npnt
        write(3,111)' Q = ',k,(pnt(i,k),i=1,3),sym_label(k)
        write(4,111)' Q = ',k,(pnt(i,k),i=1,3),sym_label(k)
        write(5,111)' Q = ',k,(pnt(i,k),i=1,3),sym_label(k)
        write(6,111)' Q = ',k,(pnt(i,k),i=1,3),sym_label(k)
        do i_nu=0,n_nu
          write(3,'(f10.4,3(1x,e15.8))')w_nu(i_nu)*fc,
     &      (w_output_int(i_nu,i,k)*fc,i=1,3)
          write(4,'(f10.4,3(1x,e15.8))')w_nu(i_nu)*fc,
     &      (p_output_int(i_nu,i,k),i=1,3)
          write(5,'(f10.4,4(1x,e15.8))')w_nu(i_nu)*fc,
     &      (e_output_int(i_nu,i,k),i=1,4)
          write(6,'(f10.4,4(1x,e15.8))')w_nu(i_nu)*fc,
     &      (e1_output_int(i_nu,i,k),i=1,3)
        enddo
      enddo
      close(3)
      close(4)
      close(5)
      close(6)
c     ------ Plot along line ---------------------------------------
      h_r1=(r1_end-r1_begin)/ndiv_r1
      hd=par*sqrt(dot_product(h_r1,h_r1))
      open(3,file=allfile(1:i_len)//'_w_nu_line')
      do idiv=0,ndiv_r1
        write(3,'(f9.5,4(1x,e15.8))')idiv*hd,
     &    (w_pr_nu(idiv,i_nu),i_nu=0,3)
      enddo
      close(3)
      open(3,file=allfile(1:i_len)//'_w_tau_line')
      do idiv=0,ndiv_r1
        write(3,'(f9.5,4(1x,e15.8))')idiv*hd,
     &    (w_pr_tau(idiv,i_tau),i_tau=0,3)
      enddo
      close(3)
c     ------- G=G'=0 components for q-points of high symmetry -----------
      allocate(bos(0:n_nu,n_symb))
      bos=0.d0
      do k=1,n_symb
        lb=ksym_in_npnt(k)
        if(lb/=0) bos(:,k)=p_output_int(:,1,lb)
      enddo
      call output_boson_q_nu(bos,'_P  ')
      bos=0.d0
      do k=1,n_symb
        lb=ksym_in_npnt(k)
        if(lb/=0) bos(:,k)=w_output_int(:,1,lb)
      enddo
      call output_boson_q_nu(bos,'_W  ')
      bos=0.d0
      do k=1,n_symb
        lb=ksym_in_npnt(k)
        if(lb/=0) bos(:,k)=e_output_int(:,1,lb)
      enddo
      call output_boson_q_nu(bos,'_E  ')
      bos=0.d0
      do k=1,n_symb
        lb=ksym_in_npnt(k)
        if(lb/=0) bos(:,k)=e_output_int(:,2,lb)
      enddo
      call output_boson_q_nu(bos,'_E_M')
      deallocate(bos)
c$$$      call e_real_axis(e_output_int)
 111  format(a5,i4,1x,3f8.5,4x,a1)      
      end
