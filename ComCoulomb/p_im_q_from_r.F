      subroutine p_im_q_from_r(p_mi,it,ind_tau,key,key_p)
c     key = 0 : contribution to P_W_SOLID_Q_TAU -----------------
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      include "mpif.h"
      integer, intent(in) :: ind_tau,key,it,key_p
      double precision, intent(in) ::
     $  p_mi(nr_red,n_pbmt,ndim3_kk(me3_k+1))
      integer :: ind,i,i0,j,ia,ib,ind_k,k,ig,jatom,jsort,
     $  ind_ir,ir,k0,k1,ii,jj
      double precision :: v(3),fc,pi2
      complex*16 :: cc, p_c2
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),tmp2(:,:),s(:,:),
     &  p_c(:,:),pcc(:,:,:,:),tmp3(:,:),
     &  tmp4(:,:)
      pi2=pi+pi
      
!      write(100,*)'ind_tau',ind_tau,'key',key,
!     & 'it',it,'key_p',key_p
!      if(maswrk) write(100,*)'gbs',gbs,'ind_k_star_back',
!     $ ind_k_star_back,'indpw_gw',indpw_gw,'indpw_r',indpw_r,
!     $ 'iopb',iopb,'ir_ref',ir_ref,'is',is,'k_list',k_list,
!     $ 'k_star',k_star,'k_star_max',k_star_max,'k_sym_0',k_sym_0,
!     $ 'n_pbmt',n_pbmt,'n_pbmt0',n_pbmt0,'n_pbtot',n_pbtot,
!     $ 'natom',natom,'ndiv',ndiv,'nplw_gw',nplw_gw,'nplw_r',nplw_r,
!     $ 'nplw_rm',nplw_rm,'nplwgw',nplwgw,'nqdiv',nqdiv,'nr_full',
!     $ nr_full,'nr_red',nr_red,'nrdiv',nrdiv,'nvopt',nvopt,'nvopt_max',
!     $ nvopt_max,'p_local_r',p_local_r,'p_w_solid_q_tau',p_w_solid_q_tau
!     $ ,'pnt',pnt,'qvopt',qvopt,'r0b0',r0b0,'r_group',r_group,'r_pnt',
!     $ r_pnt,'rb0',rb0,'rbas',rbas,'tshift',tshift,'tshift_r',tshift_r,
!     $ 'v_repr',v_repr,'xi_0_q_tau',xi_0_q_tau
      
      allocate(p_c(nqdiv,n_pbmt))
      allocate(pcc(n_pbmt,nr_red,k_star_max,ndim3_k(me3_k+1)))
      
      do i0=1,nr_red
        p_c=(0.d0,0.d0)
        do ind_ir=1,ndim3_kk(me3_k+1)
          ir=n3_mpi_kk(me3_k+1)+ind_ir
          do j=1,n_pbmt
            if ((key_p .ne. 0) .and. (p_local_r)) then
              if (ir .eq. 1) then
                p_c(ir,j)=p_mi(i0,j,ind_ir)
              endif
            else
              p_c(ir,j)=p_mi(i0,j,ind_ir)
            endif
          enddo
        enddo
        if(nproc_k/=1) call mpi_allreduce_dcmplx
     $    (p_c,nqdiv*n_pbmt,mpi_sum,comm_tau)
        
c$$$        if ((key_p .eq. 2) .and. (myrevision .eq. 1)) then
c$$$          do i=1,n_pbmt
c$$$            p_c2=p_c(1,i)
c$$$            do k=1,nqdiv
c$$$              p_c(k,i)=p_c2
c$$$            enddo
c$$$          enddo
c$$$        else
!        if(maswrk) write(100,*)'bk2'
        call fft3(ndiv(1),ndiv(2),ndiv(3),n_pbmt,p_c,1)
c$$$        endif
!        if(maswrk) write(100,*)'bk3'
        do ind_k=1,ndim3_k(me3_k+1)
          k=n3_mpi_k(me3_k+1)+ind_k
          do k0=1,k_star(k)
            k1=k_list(k0,k)
            call zone1_number(pnt(1,k1),rb0,ndiv,ind)
            pcc(:,i0,k0,ind_k)=p_c(ind,:)
          enddo
        enddo
      enddo                     !! over i0
      deallocate(p_c)
!      if(maswrk) write(100,*)'bk3'
      allocate(tmp(nr_full,n_pbmt))
      allocate(tmp1(n_pbmt,nplw_rm))
      if(key==0) then
        allocate(tmp3(nvopt_max,nplw_rm))
        allocate(tmp4(nvopt_max,nvopt_max))
      endif
!      if(v_repr=='NON')	allocate(tmp2(n_pbmt,nplw_gw))
      allocate(tmp2(n_pbmt,nplw_gw))
!      if(maswrk) write(100,*)'bk4'
c     ---------------------------------------------------------------	
      do ind_k=1,ndim3_k(me3_k+1)
        k=n3_mpi_k(me3_k+1)+ind_k
c     ------- We fill in the first index completely -------------
        do i=1,nr_full
          i0=ir_ref(i)
          ig=r_group(i)
          ind=k_sym_0(k,ig)     ! IND=A^(-1)*K  in NQDIV-list
          ind=ind_k_star_back(ind) ! position in a given K-star
          call sym_vmi_left(ig,pcc(:,i0,ind,ind_k),tmp(i,:))
          do jatom=1,natom
            jsort=is(jatom)
            v=tshift(:,jatom,ig)-tshift_r(:,i,ig)
            fc=pi2*dot_product(pnt(:,k),v)
            cc=dcmplx(dcos(fc),dsin(fc))
            ia=iopb(jatom)
            ib=ia+n_pbmt0(jsort)-1
            tmp(i,ia:ib)=tmp(i,ia:ib)*cc
          enddo
c     ----- Multiplication with the phase factors ------------------
          v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &      +r_pnt(3,i)*r0b0(:,3)
          fc=2.d0*pi*dot_product(pnt(:,k),v)
          cc=cdexp(dcmplx(0.d0,fc))
          tmp(i,1:n_pbmt)=tmp(i,1:n_pbmt)*cc
        enddo                   !! over ist
c     -----------------------------------------------------------
        call fft3(nrdiv(1),nrdiv(2),nrdiv(3),n_pbmt,tmp,1)
        tmp=tmp/nr_full
        tmp1=(0.d0,0.d0)
          do i=1,nplw_r(k)
            i0=indpw_r(i,k)
            call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
            tmp1(1:n_pbmt,i)=tmp(ind,1:n_pbmt)
          enddo
c     ------ We form the overlap matrix -----------------------------
          allocate(s(nplw_r(k),nplwgw(k)))
          call s_overlap_m_b(s,nplw_r(k),nplwgw(k),indpw_r(1,k),
     &      indpw_gw(1,k))
          call zgemm('n','n',n_pbmt,nplwgw(k),nplw_r(k),(1.d0,0.d0),
     &      tmp1,n_pbmt,s,nplw_r(k),(0.d0,0.d0),tmp2,n_pbmt)
          tmp1(:,1:nplwgw(k))=tmp2(:,1:nplwgw(k))
          deallocate(s)

c$$$        do j=1,nplwgw(k)
c$$$          do i=1,n_pbmt
c$$$            p_w_solid_q_tau2(i,n_pbmt+j,ind_tau,ind_k)
c$$$     $        =p_w_solid_q_tau2(i,n_pbmt+j,ind_tau,ind_k)
c$$$     &        +dreal(tmp1(i,j))
c$$$            p_w_solid_q_tau2(n_pbmt+j,i,ind_tau,ind_k)
c$$$     $        =p_w_solid_q_tau2(n_pbmt+j,i,ind_tau,ind_k)
c$$$     &        +dimag(tmp1(i,j))
c$$$          enddo
c$$$        enddo
!        if(maswrk) write(100,*)'bk5'
        if(key==0) then

c$$$          write(179,*) 'p_w_solid_q_tau2, im'          
c$$$          do ii=1, n_pbmt
c$$$            do jj=1, nplw_rm
c$$$              p_w_solid_q_tau2(ii,n_pbmt+jj,ind_tau,ind_k,key_p)
c$$$     $          =tmp1(ii,jj)                    
c$$$              p_w_solid_q_tau2(n_pbmt+jj,ii,ind_tau,ind_k,key_p)
c$$$     $          =dconjg(tmp1(ii,jj))
c$$$            enddo
c$$$          enddo
          call zgemm('c','n',nvopt(k),nplwgw(k),n_pbmt,(1.d0,0.d0),
     &      qvopt(1,1,ind_k),n_pbtot,tmp1,n_pbmt,(0.d0,0.d0),
     &      tmp3,nvopt_max)
!        if(maswrk) write(100,*)'bk6'
          call zgemm('n','n',nvopt(k),nvopt(k),nplwgw(k),(1.d0,0.d0),
     &      tmp3,nvopt_max,qvopt(n_pbmt+1,1,ind_k),n_pbtot,
     &      (0.d0,0.d0),tmp4,nvopt_max)
!        if(maswrk) write(100,*)'bk7'
          do j=1,nvopt(k)
            do i=1,j
              tmp4(i,j)=tmp4(i,j)+conjg(tmp4(j,i))
              if(i/=j) tmp4(j,i)=conjg(tmp4(i,j))
            enddo
          enddo
          call boson_pack_tau(1.d0,tmp4,nvopt_max,1.d0,
     &      p_w_solid_q_tau(1,1,ind_tau,ind_k),
     &      nvopt_max,nvopt(k))
        else if(key==1) then
          do j=1,nplwgw(k)
            do i=1,n_pbmt
              xi_0_q_tau(i,n_pbmt+j,it,ind_tau,ind_k)=dble(tmp1(i,j))
              xi_0_q_tau(n_pbmt+j,i,it,ind_tau,ind_k)=dimag(tmp1(i,j))
            enddo
          enddo
        endif
c     --------------------------------------------------------------
!        if(maswrk) write(100,*)'bk8'
      enddo                     !! over ind_k
!      if(v_repr=='NON')	deallocate(tmp2)
      deallocate(tmp2)
      deallocate(tmp,tmp1,pcc)
      if(key==0) deallocate(tmp3,tmp4)
      end
