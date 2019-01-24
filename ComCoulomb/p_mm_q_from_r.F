      subroutine p_mm_q_from_r(p_mm,it,ind_tau,key,key_p)

c     key = 0 : contribution to P_W_SOLID_Q_TAU -----------------
c     key = 1 : contribution to XI_0_Q_TAU -----------------

      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod

      implicit none
      include "mpif.h"
      integer, intent(in) :: ind_tau,key,it,key_p
      double precision, intent(in) ::
     $  p_mm(n_pbmt,n_pbmt,ndim3_kk(me3_k+1))
      integer :: ind,ind_k,k,i,j,ind_ir,ir,ii,jj
      double precision :: dev0
      complex*16 :: p_c2
      double precision, allocatable :: tmp(:,:,:)
      complex*16, allocatable :: p_c(:,:),tmp1(:,:),tmp2(:,:),tmp3(:,:)


      if(key==0) allocate(tmp(n_pbmt,n_pbmt,ndim3_k(me3_k+1)))
      allocate(p_c(nqdiv,n_pbmt))



      do j=1,n_pbmt
        p_c=(0.d0,0.d0)
        do ind_ir=1,ndim3_kk(me3_k+1)
          ir=n3_mpi_kk(me3_k+1)+ind_ir
          do i=1,j
            if ((key_p .ne. 0) .and. (p_local_r)) then
              if (ir .eq. 1) then
                p_c(ir,i)=p_mm(i,j,ind_ir)
              endif
            else
              p_c(ir,i)=p_mm(i,j,ind_ir)
            endif
          enddo
        enddo
        
        if(nproc_k/=1) call mpi_allreduce_dcmplx
     $    (p_c,nqdiv*j,mpi_sum,comm_tau)
c$$$  if ((key_p .eq. 2) .and. (myrevision .eq. 1)) then
c$$$  do i=1,j
c$$$  p_c2=p_c(1,i)
c$$$  do k=1,nqdiv
c$$$  p_c(k,i)=p_c2
c$$$  enddo
c$$$  enddo
c$$$  else
        call fft3(ndiv(1),ndiv(2),ndiv(3),j,p_c,1)          
c$$$  endif
        do ind_k=1,ndim3_k(me3_k+1)
          k=n3_mpi_k(me3_k+1)+ind_k
          call zone1_number(pnt(:,k),rb0,ndiv,ind)
c$$$          do i=1,j-1
c$$$	      p_w_solid_q_tau2(i,j,ind_tau,ind_k)
c$$$     $        =p_w_solid_q_tau2(i,j,ind_tau,ind_k)+real(p_c(ind,i))
c$$$	      p_w_solid_q_tau2(j,i,ind_tau,ind_k)
c$$$     $          =p_w_solid_q_tau2(j,i,ind_tau,ind_k)+imag(p_c(ind,i))
c$$$          enddo
c$$$          p_w_solid_q_tau2(j,j,ind_tau,ind_k)
c$$$     $      =p_w_solid_q_tau2(j,j,ind_tau,ind_k)+real(p_c(ind,j))
          if(key==0) then
            do i=1,j-1
              tmp(i,j,ind_k)=2.d0*dble(p_c(ind,i))
              tmp(j,i,ind_k)=2.d0*dimag(p_c(ind,i))
            enddo
            tmp(j,j,ind_k)=2.d0*dble(p_c(ind,j))
          else if(key==1) then
            do i=1,j-1
              xi_0_q_tau(i,j,it,ind_tau,ind_k)=dble(p_c(ind,i))
              xi_0_q_tau(j,i,it,ind_tau,ind_k)=dimag(p_c(ind,i))
            enddo
            xi_0_q_tau(j,j,it,ind_tau,ind_k)=dble(p_c(ind,j))
          endif
        enddo                   !! over ind_k
      enddo                     !! over j
        
      deallocate(p_c)
      
      if(key==0) then
        allocate(tmp1(n_pbmt,n_pbmt))
        allocate(tmp2(nvopt_max,n_pbmt))
        allocate(tmp3(nvopt_max,nvopt_max))
        do ind_k=1,ndim3_k(me3_k+1)
          k=n3_mpi_k(me3_k+1)+ind_k
          tmp1=0.0d0
          call boson_unpack_tau(tmp1,n_pbmt,tmp(1,1,ind_k),n_pbmt,
     &      n_pbmt)
c$$$          write(179,*) 'p_w_solid_q_tau2, mm'
c$$$          do ii=1, n_pbmt
c$$$            do jj=1, n_pbmt
c$$$              p_w_solid_q_tau2(ii,jj,ind_tau,ind_k,key_p)=tmp1(ii,jj)
c$$$            enddo
c$$$          enddo
          call zgemm('c','n',nvopt(k),n_pbmt,n_pbmt,(1.d0,0.d0),
     &      qvopt(1,1,ind_k),n_pbtot,tmp1,n_pbmt,(0.d0,0.d0),
     &      tmp2,nvopt_max)
          call zgemm('n','n',nvopt(k),nvopt(k),n_pbmt,(1.d0,0.d0),
     &      tmp2,nvopt_max,qvopt(1,1,ind_k),n_pbtot,
     &      (0.d0,0.d0),tmp3,nvopt_max)
          call boson_pack_tau(1.d0,tmp3,nvopt_max,1.d0,
     &      p_w_solid_q_tau(1,1,ind_tau,ind_k),
     &      nvopt_max,nvopt(k))
        enddo                   !! over ind_k

        deallocate(tmp,tmp1,tmp2,tmp3)
      endif
      end
