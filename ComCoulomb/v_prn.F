      subroutine v_prn
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      include "mpif.h"
      integer :: iq,iq0,ndim0,j,i,i_len,idiv,i_start,ndim,iat_r0,ind_k,
     &  n_pbt,kst,key
      double precision :: r_transl(3),h_r1(3),hd,d0,r0(3)
      character*2 :: ubi_r
      complex*16 :: zdotc
      integer, allocatable :: iat_r(:)
      double precision, allocatable :: v_pr(:)
      complex*16, allocatable :: w_0(:,:),work(:),y(:),v_tmp(:,:),
     &  s1(:,:)
!      if(v_repr=='NON') key=0
      key=0
      allocate(w_0(n_pbtot,n_pbtot))
      allocate(y(n_pbtot))
      allocate(work(n_pbtot))
      allocate(iat_r(0:ndiv_r1))
      allocate(v_pr(0:ndiv_r1))
      h_r1=(r1_end-r1_begin)/ndiv_r1
      ndim0=max(n_pbmtm,nplw_gw)
      v_pr=0.d0
      do ind_k=1,ndim3_k(me3_k+1)
        iq0=n3_mpi_k(me3_k+1)+ind_k
        n_pbt=n_pbmt+nplwgw(iq0)
        allocate(v_tmp(n_pbt,n_pbt))
        call v_coul_full(nplwgw(iq0),v_tmp,ind_k) 
        do kst=1,k_star(iq0)
          iq=k_list(kst,iq0)
          do j=1,n_pbt
            do i=1,n_pbt
              w_0(i,j)=v_tmp(i,j)
            enddo
          enddo
c     -------- Obtaining V in the point Q from V in Q0 -----------------
          call sym_w(iq,iq0,w_0)
c     ----------- Product basis at the point R0 -----------------------
          allocate(s1(nplwgw(iq0),nplwgw(iq0)))
          if(key==1) then
            if(nplwgw(iq0)/=0) call s_overlap(1,s1,nplwgw(iq0),
     &        nplwgw(iq0),
     &        indpw_gw(1,iq))
          endif
          call value_pb(r_ref,iq,iq0,ubi_r,iat_r0,r_transl,work,i_start,
     &      ndim,ndim0,key,s1,nplwgw(iq0))
c$$$          write(179,*) 'ind_k, ndim', ndim, n_pbtot
c$$$          call flush(179)
          if(ndim/=0) call zgemv('t',ndim,n_pbt,(1.d0,0.d0),
     &      w_0(i_start,1),n_pbtot,work,1,
     &      (0.d0,0.d0),y,1)
          do idiv=0,ndiv_r1
            call value_pb(r1_begin+idiv*h_r1,iq,iq0,ubi_r,iat_r(idiv),
     &        r_transl,work,i_start,ndim,ndim0,key,s1,
     &        nplwgw(iq0))
            if(ndim/=0) v_pr(idiv)=v_pr(idiv)
     &        +zdotc(ndim,work,1,y(i_start),1)/nqdiv
          enddo
          deallocate(s1)
        enddo                   !! over iq
        deallocate(v_tmp)
      enddo                     !! over ind_k
      deallocate(w_0,y,work)
      if(nproc_k/=1) call mpi_allreduce_dble
     $  (v_pr,ndiv_r1+1,mpi_sum,comm_tau)
      hd=par*sqrt(dot_product(h_r1,h_r1))
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_v_coul_line')
        do idiv=0,ndiv_r1
          r0=r1_begin+idiv*h_r1-r_ref
          d0=par*sqrt(dot_product(r0,r0))
          if(d0<1.d-10) cycle
          write(3,'(2i3,4(e15.8,1x))')iat_r0,iat_r(idiv),idiv*hd,
     &      v_pr(idiv),2.d0/d0,v_pr(idiv)-2.d0/d0
        enddo
        close(3)
      endif
      deallocate(v_pr,iat_r)
      end
