      subroutine p_ii_solid(it0,ind_tau,key,key1,key_p)

c     key = 0 : (G_tau)*(G_beta-tau) contribution to P -----------------
c     key = 1 : (G_tau_UP)*(G_beta-tau_DN) contribution to Xi_T -----
c     key1 = 0 - to use QP Green's function
c     key1 = 1 - to use regular Green's function

      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod


      implicit none
      include "mpif.h"
      integer, intent(in) :: ind_tau,key,it0,key1,key_p
      integer :: ind_s,nsta,ispin,ind_j,ind_i,i,ir,jrl,irl,ind_k,k,i0,
     &  ig,ind,j,iv(3),j0,j0s,k0,k1,it,nsp,ii,jj
      double precision :: v(3),fc,pi2,dum,dev0
      complex*16 :: cc,dumc,p_c2
      double precision, allocatable :: gr1(:,:,:),gr2(:,:,:)
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),p_c(:,:,:),
     &  p_tmp(:,:,:),g_c(:,:,:,:,:,:),
     &  gc1(:,:,:,:,:),gc2(:,:,:,:,:),s(:,:),
     &  pcc(:,:,:,:),pc(:,:),tmp2(:,:),
     &  tmp3(:,:), tmp4(:,:)

!      write(100,*)'gbs',gbs,'ind_k_star_back',ind_k_star_back,
!     $'ind_r_star',ind_r_star,'ind_r_star_back',ind_r_star_back,
!     $'indpw_gw',indpw_gw,'indpw_r',indpw_r,'ip_r',ip_r,'ir_ref',
!     $ir_ref,'irel',irel,'itshift_r',itshift_r,'k_list',k_list,
!     $'k_star',k_star,'k_star_max',k_star_max,'k_sym_0',k_sym_0,
!     $'max_r_star',max_r_star,'n_pbmt',n_pbmt,'n_pbtot',n_pbtot,
!     $'nbasmpw',nbasmpw,'ndiv',ndiv,'nplw_r',nplw_r,'nplw_rm',
!     $nplw_rm,'nplwgw',nplwgw,'nqdiv',nqdiv,'nr_full',nr_full,
!     $'nr_red',nr_red,'nrdiv',nrdiv,'nrel',nrel,'nspin',nspin,
!     $'num_r_star',num_r_star,'nvopt',nvopt,'nvopt_max',nvopt_max,
!     $'p_local_r',p_local_r,'p_w_solid_q_tau',p_w_solid_q_tau,
!     $'phase_k_rr',phase_k_rr,'pnt',pnt,'qvopt',qvopt,'r0b0',r0b0,
!     $'r_group',r_group,'r_pnt',r_pnt,'rb0',rb0,'rbas',rbas,
!     $'v_repr',v_repr,'xi_0_q_tau',xi_0_q_tau 

!      call mpi_bcast(natom,1,mpi_integer,0,mpi_comm_world,mpierr)

      if(key==0) nsp=nspin
      if(key==1) nsp=1
      pi2=pi+pi
      allocate(g_c(nbasmpw,ndim3_k(me3_k+1),nr_full,nrel,2,nsp))
      if(key==0) then
        do it=1,2
          call g_k_g_r1(it,ind_tau,g_c,key1,key_p)
        enddo
      else if(key==1) then
        call g_k_g_r1_xi(it0,ind_tau,g_c(1,1,1,1,1,1),1,key1)
        call g_k_g_r1_xi(3-it0,ind_tau,g_c(1,1,1,1,2,1),nspin,key1)
      endif

c      call log_eachcore('** g_k_g_r1 done        **********')              
      allocate(p_tmp(nplw_rm,ndim3_k(me3_k+1),nr_full))
      allocate(p_c(nqdiv,ndim3_nr_red(me3_k+1),max_r_star))
      allocate(pcc(nr_red,max_r_star,k_star_max,ndim3_k(me3_k+1)))		
      allocate(tmp(nr_full,max_r_star))
      allocate(pc(nqdiv,nr_red))
      if(irel<=1) then
        allocate(gr1(nqdiv,ndim3_nr_red(me3_k+1),max_r_star))
        allocate(gr2(nqdiv,ndim3_nr_red(me3_k+1),max_r_star))
      else if(irel==2) then
        allocate(gc1(nqdiv,ndim3_nr_red(me3_k+1),max_r_star,nrel,nrel))
        allocate(gc2(nqdiv,ndim3_nr_red(me3_k+1),max_r_star,nrel,nrel))
      endif
      p_tmp=(0.d0,0.d0)
c     --- P_C (R,r,r') where R - all; r - irred; r' - current star
c$$$  call log_eachcore('** loop start        **********')                    
      do ind_s=1,nr_red         !! over stars of r'
        nsta=num_r_star(ind_s)  !! size of the star
        p_c=(0.d0,0.d0)
        do ispin=1,nsp
          if(irel<=1) then
            call g_rr_r_r1_star(ind_s,nsta,g_c(1,1,1,1,1,ispin),gr1,
     &        dumc)
            call g_rr_r_r1_star(ind_s,nsta,g_c(1,1,1,1,2,ispin),gr2,
     &        dumc)
          else if(irel==2) then
            call g_rr_r_r1_star(ind_s,nsta,g_c(1,1,1,1,1,ispin),dum,gc1)
            call g_rr_r_r1_star(ind_s,nsta,g_c(1,1,1,1,2,ispin),dum,gc2)
          endif  
          do ind_j=1,nsta       !! over r' star IND_S
            do ind_i=1,ndim3_nr_red(me3_k+1) !! irreducible r
              if(irel<=1) then
                do ir=1,nqdiv
                  p_c(ir,ind_i,ind_j)=p_c(ir,ind_i,ind_j)
     &              -gr1(ir,ind_i,ind_j)
     &              *gr2(ir,ind_i,ind_j)
                enddo
              else if(irel==2) then
                do jrl=1,2
                  do irl=1,2
                    do ir=1,nqdiv
                      p_c(ir,ind_i,ind_j)=p_c(ir,ind_i,ind_j)
     &                  -gc1(ir,ind_i,ind_j,irl,jrl)
     &                  *conjg(gc2(ir,ind_i,ind_j,irl,jrl))
                    enddo
                  enddo
                enddo           !! over jrl
              endif
            enddo               !! over ind_i
          enddo                 !! over ind_j
        enddo                   !! over ispin

        if ((key_p .ne. 0) .and. (p_local_r)) then
          do ir=2,nqdiv
            do ind_j=1,nsta     !! over r' star IND_S
              do ind_i=1,ndim3_nr_red(me3_k+1) !! irreducible r
                p_c(ir,ind_i,ind_j)=0.0d0
              enddo
            enddo
          enddo
        endif        
        

        
        do ind_j=1,nsta         !! over r' star IND_S
          pc=(0.d0,0.d0)
          do ind_i=1,ndim3_nr_red(me3_k+1) !! irreducible r
            i=n3_mpi_nr_red(me3_k+1)+ind_i
c$$$            if ((key_p .eq. 2) .and. (myrevision .eq. 1)) then
c$$$              p_c2=p_c(1,ind_i,ind_j)
c$$$              do ii=1, nqdiv
c$$$                p_c(ii,ind_i,ind_j)=p_c2
c$$$              enddo
c$$$            else
            call fft3(ndiv(1),ndiv(2),ndiv(3),1,p_c(1,ind_i,ind_j),1)
c$$$            endif
            pc(:,i)=p_c(:,ind_i,ind_j)
          enddo
          if(nproc_k/=1) call mpi_allreduce_dcmplx
     $      (pc,nqdiv*nr_red,mpi_sum,comm_tau)
          do ind_k=1,ndim3_k(me3_k+1)
            k=n3_mpi_k(me3_k+1)+ind_k
            do k0=1,k_star(k)
              k1=k_list(k0,k)
              call zone1_number(pnt(1,k1),rb0,ndiv,ind)
              pcc(:,ind_j,k0,ind_k)=pc(ind,:)
            enddo
          enddo
        enddo

c$$$        call log_eachcore('** loop end       **********')                            
c     ---------------------------------------------------------------
        do ind_k=1,ndim3_k(me3_k+1)
          k=n3_mpi_k(me3_k+1)+ind_k
          tmp(:,1:nsta)=(0.d0,0.d0)
c     ----------- Restore all r --------------------------------          
          do i=1,nr_full
            i0=ir_ref(i)
            ig=r_group(i)
            ind=k_sym_0(k,ig)   ! IND=A^(-1)*K  in NQDIV-list
            ind=ind_k_star_back(ind) ! position in a given K-star
            do ind_j=1,nsta
              j=ind_r_star(ind_j,ind_s)
              iv=itshift_r(:,i,ig)-itshift_r(:,j,ig)
              j0=ip_r(j,ig)
              j0s=ind_r_star_back(j0)
              tmp(i,ind_j)=pcc(i0,j0s,ind,ind_k)
     &          *phase_k_rr(k,iv(1),iv(2),iv(3))
            enddo
          enddo
c     -----------------------------------------------------------	    
          do i=1,nr_full
            v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &        +r_pnt(3,i)*r0b0(:,3)
            fc=pi2*dot_product(pnt(:,k),v)
            cc=exp(dcmplx(0.d0,-fc))
            tmp(i,1:nsta)=tmp(i,1:nsta)*cc
          enddo
          call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nsta,tmp,-1)
          do ind_j=1,nsta
            tmp(:,ind_j)=tmp(:,ind_j)/nr_full
          enddo
c     ------------------------------------------------------------
            do i=1,nplw_r(k)
              i0=indpw_r(i,k)
              call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
              do ind_j=1,nsta
                j=ind_r_star(ind_j,ind_s)
                p_tmp(i,ind_k,j)=tmp(ind,ind_j)
              enddo
            enddo
        enddo                   !! over ind_k
      enddo                     !! over ind_s

c$$$      call log_eachcore('** loop 2nd end       **********')                                  
      deallocate(p_c,tmp,pcc,pc)
      if(irel<=1) deallocate(gr1,gr2)
      if(irel==2) deallocate(gc1,gc2)
      allocate(tmp(nr_full,nplw_rm))
      allocate(tmp1(nplw_rm,nplw_rm))
c$$$      allocate(tmp4(nplw_rm,nplw_rm))      
      if(key==0) then
        allocate(tmp2(nvopt_max,nplw_rm))
        allocate(tmp3(nvopt_max,nvopt_max))
      endif
      do ind_k=1,ndim3_k(me3_k+1)
        k=n3_mpi_k(me3_k+1)+ind_k
          do i=1,nr_full
            v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &        +r_pnt(3,i)*r0b0(:,3)
            fc=pi2*dot_product(pnt(:,k),v)
            cc=exp(dcmplx(0.d0,fc))
            tmp(i,1:nplw_r(k))=p_tmp(1:nplw_r(k),ind_k,i)*cc
          enddo
          call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nplw_r(k),tmp,1)
          do i=1,nplw_r(k)
            i0=indpw_r(i,k)
            call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
            do j=1,nplw_r(k)
              tmp1(j,i)=tmp(ind,j)/nr_full
            enddo
          enddo
c     ------ We form the overlap matrix -----------------------------
          allocate(s(nplw_r(k),nplwgw(k)))
          call s_overlap_m_b(s,nplw_r(k),nplwgw(k),indpw_r(1,k),
     &      indpw_gw(1,k))
          call zgemm('n','n',nplw_r(k),nplwgw(k),nplw_r(k),
     &      (1.d0,0.d0),tmp1,nplw_rm,s,nplw_r(k),(0.d0,0.d0),
     &      tmp,nr_full)
          call zgemm('C','n',nplwgw(k),nplwgw(k),nplw_r(k),
     &      (1.d0,0.d0),s,nplw_r(k),tmp,nr_full,(0.d0,0.d0),
     &      tmp1,nplw_rm)
          deallocate(s)

c$$$        do j=1,nplwgw(k)
c$$$          do i=1,j
c$$$            p_w_solid_q_tau2(n_pbmt+i,n_pbmt+j,ind_tau,ind_k)
c$$$     $        =p_w_solid_q_tau2(n_pbmt+i,n_pbmt+j,ind_tau,ind_k)
c$$$     $        +dreal(tmp1(i,j))
c$$$            if(i/=j) p_w_solid_q_tau2(n_pbmt+j,n_pbmt+i,ind_tau,ind_k)
c$$$     $        =p_w_solid_q_tau2(n_pbmt+j,n_pbmt+i,ind_tau,ind_k)
c$$$     $        +dimag(tmp1(i,j))
c$$$          enddo
c$$$        enddo
        
        if(key==0) then
c$$$          write(179,*) 'p_w_solid_q_tau2, ii'          
c$$$          do ii=1, nplw_rm
c$$$            do jj=1, nplw_rm            
c$$$              p_w_solid_q_tau2
c$$$     $          (n_pbmt+ii,n_pbmt+jj,ind_tau,ind_k,key_p)
c$$$     $          =tmp1(ii,jj)
c$$$            enddo
c$$$          enddo
          call zgemm('c','n',nvopt(k),nplwgw(k),nplwgw(k),(1.d0,0.d0),
     &      qvopt(n_pbmt+1,1,ind_k),n_pbtot,tmp1,nplw_rm,
     &      (0.d0,0.d0),tmp2,nvopt_max)
          call zgemm('n','n',nvopt(k),nvopt(k),nplwgw(k),(1.d0,0.d0),
     &      tmp2,nvopt_max,qvopt(n_pbmt+1,1,ind_k),n_pbtot,
     &      (0.d0,0.d0),tmp3,nvopt_max)
          call boson_pack_tau(1.d0,tmp3,nvopt_max,1.d0,
     &      p_w_solid_q_tau(1,1,ind_tau,ind_k),
     &      nvopt_max,nvopt(k))
        else if(key==1) then
          do j=1,nplwgw(k)
            do i=1,j-1
              xi_0_q_tau(n_pbmt+i,n_pbmt+j,it0,ind_tau,ind_k)=
     &          dble(tmp1(i,j))
              xi_0_q_tau(n_pbmt+j,n_pbmt+i,it0,ind_tau,ind_k)=
     &          dimag(tmp1(i,j))
            enddo
            xi_0_q_tau(n_pbmt+j,n_pbmt+j,it0,ind_tau,ind_k)=
     &        dble(tmp1(j,j))
          enddo
        endif
      enddo                     !! over ind_k
c      call log_eachcore('** the rest end      **********')                                        
      deallocate(tmp1,tmp,p_tmp)
      if(key==0) deallocate(tmp2,tmp3)
      end
