      subroutine g_k_g_r1(it,ind_tau,g_c,key,key_p)

c$$$  Calculates G^k(G,r_1) where k=1,npnt+par; G=1,nbasmpw; r'=1,nr_full -
c$$$  where k is in the IBZ
c$$$  key = 0 - to use QP Green's function
c$$$  key = 1 - to use regular Green's function
      
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use wannier90_mod
      
      implicit none
      integer, intent(in) :: ind_tau,it,key,key_p
      complex*16, intent(out) :: g_c(nbasmpw,ndim3_k(me3_k+1),nr_full,
     &  nrel,2,nspin)
      integer :: k,nbnd,ind_k,nbas,irl,ibas0,i,ispin,i0,ind,i_tau,j,ib,
     &  jb,ii,jj
      double precision :: v(3),fc,pi2,tu,de,g_diag,green0
      complex*16 :: cc
      double precision, allocatable :: gx_tau(:)
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),tmp0(:,:),
     &  tmp2(:,:),tempvec(:,:)

      allocate(tmp0(nbndf,nbasmpw),tmp2(nbndf,nbndf))
      allocate(tmp1(nr_full,nbasmpw))
      allocate(tmp(nbasmpw,nbasmpw))

      pi2=pi+pi
      allocate(gx_tau(nbndf))
c     ---------------------------------------------------------------------

c$$$  if (maswrk) then
c$$$  write(993, *) 'indgb(nbndtm, nqdiv)'
c$$$  do ii=1, nbndtm
c$$$  do jj=1, nqdiv
c$$$  write(993,'(3(i5, 2x), 3f12.6)') ii, jj,
c$$$  $        indgb(ii, jj), gbs(:,indgb(ii, jj))
c$$$  enddo
c$$$  enddo
c$$$  endif
      
      i_tau=me3_tau*ndim3_tau+ind_tau-1
      if(it==2) i_tau=n_tau-i_tau      
      
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me3_k+1)
          k=n3_mpi_k(me3_k+1)+ind_k
          nbnd=n_bnd(k,ispin)

c$$$  start1 : get \sum_{\lamdba,\lambda'}A_{G}^{k,\lamdba} G_{\lambda,lambda'}^{k}(\tau) A_{G'}^{* k,\lamdba'} in eq. 22 in IBZ
          
          if(ubi=='dft'.or.ubi=='dif'.or.ubi==' hf') then
            tu=tau_mesh(i_tau)
            call g_x_tau(ispin,k,gx_tau,tu,nbnd,chem_pot)
            do j=1,nbask(k)
              do ib=1,nbnd
                tmp0(ib,j)=gx_tau(ib)*conjg(ev_bnd(j,ib,ind_k,ispin))
              enddo
            enddo
          else
            if(key==0) then
              if (key_p .eq. 0) then
                tmp2=(0.d0,0.d0)
c$$$  basis rotation from n(energy eigenstate) to band basis
                do i=1,nbnd
                  de=e_qp(i,k,ispin)-chem_pot_qp
                  g_diag=green0(de,tau_mesh(i_tau))
                  do jb=1,nbnd
                    do ib=1,nbnd
                      tmp2(ib,jb)=tmp2(ib,jb)+q_qp(ib,i,ind_k,ispin)
     $                  *g_diag*conjg(q_qp(jb,i,ind_k,ispin))
                    enddo
                  enddo
                enddo
              else
                if (u_mode .eq. 'bnse') then
                  tmp2=(0.d0,0.d0)
c$$$  basis rotation from n(energy eigenstate) to band basis
                  do i=1, nbnd
                    if (include_bands_foru
     $                (i,ind_k,ispin,it,key_p) .eq. 1) then                
                      de=e_qp(i,k,ispin)-chem_pot_qp
                      g_diag=green0(de,tau_mesh(i_tau))
                      do jb=1,nbnd
                        do ib=1,nbnd
                          tmp2(ib,jb)=tmp2(ib,jb)
     $                      +q_qp(ib,i,ind_k,ispin)*g_diag
     &                      *conjg(q_qp(jb,i,ind_k,ispin))
                        enddo
                      enddo
                    endif
                  enddo                
                elseif (u_mode .eq. 'wnds') then
                  tmp2=(0.d0,0.d0)
c$$$  basis rotation from n(energy eigenstate) to band basis
                  do i=1, nbnd
                    if (include_bands_foru
     $                (i,ind_k,ispin,it,key_p) .eq. 1) then                
                      de=e_qp_dis(i,k,ispin)-chem_pot_qp
                      g_diag=green0(de,tau_mesh(i_tau))
                      do jb=1,nbnd
                        do ib=1,nbnd
                          tmp2(ib,jb)=tmp2(ib,jb)
     $                      +q_qp_dis(ib,i,ind_k,ispin)*g_diag
     &                      *conjg(q_qp_dis(jb,i,ind_k,ispin))
                        enddo
                      enddo
                    endif
                  enddo                    
                  
                elseif (u_mode .eq. 'enpj') then
                  tmp2=(0.d0,0.d0)
c$$$  basis rotation from n(energy eigenstate) to band basis
                  do i=1,nbnd
                    de=e_qp(i,k,ispin)-chem_pot_qp
                    if ((de .gt. enpj_emin) .and. (de .lt. enpj_emax))
     $                then
                      g_diag=green0(de,tau_mesh(i_tau))
                      do jb=1,nbnd
                        do ib=1,nbnd
                          tmp2(ib,jb)=tmp2(ib,jb)+q_qp(ib,i,ind_k,ispin)
     $                      *g_diag*conjg(q_qp(jb,i,ind_k,ispin))
                        enddo
                      enddo
                    endif
                  enddo
                endif
              endif
c$$$  endif                
            else if(key==1) then       
              call ferm_unpack_tau(tmp2,
     &          g_full(1,1,1,ind_tau,ind_k,ispin),
     &          nbnd,nbndf,nbndf,it)
            endif
            
            call zgemm('n','c',nbnd,nbask(k),nbnd,(1.d0,0.d0),tmp2,
     &        nbndf,ev_bnd(1,1,ind_k,ispin),nbasmpw,
     &        (0.d0,0.d0),tmp0,nbndf)
          endif
          call zgemm('n','n',nbask(k),nbask(k),nbnd,(1.d0,0.d0),
     &      ev_bnd(1,1,ind_k,ispin),nbasmpw,tmp0,nbndf,
     &      (0.d0,0.d0),tmp,nbasmpw)

          nbas=nbask(k)/nrel

c$$$  end 1
          do irl=1,nrel
            ibas0=(irl-1)*nbas
c     ----- Fill in the first zone for second index -----------------------
            tmp1=(0.d0,0.d0)
            do i=1,nbas
              i0=indgb(i,k)
              call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
              tmp1(ind,1:nbask(k))=tmp(1:nbask(k),ibas0+i)
            enddo
            call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nbask(k),tmp1,-1)
c     ----- Multiplication with the phase factor ------------------
            do i=1,nr_full
              v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &          +r_pnt(3,i)*r0b0(:,3)
              fc=-pi2*dot_product(pnt(:,k),v)
              cc=cdexp(dcmplx(0.d0,fc))/amega
              g_c(1:nbask(k),ind_k,i,irl,it,ispin)=tmp1(i,1:nbask(k))
     &          *cc
            enddo
          enddo                 !! over irl
        enddo                   !! over ind_k
      enddo                     !! over ispin
      deallocate(tmp0,tmp2,tmp,tmp1,gx_tau)
      end
