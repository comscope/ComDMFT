      subroutine sigc_ii_solid(ind_tau)
c     ------ (G_tau)*(W_tau) contribution to Sig_C -----------------
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ind_tau
	integer :: ind,i,i0,j,j0,ind_k,k,ig,nbnd,ir,j0s,ind_j,iv(3),
     &	         nsta,irl,jrl,ispin,it,ista,ierr,idum,ind_i
	real*8 :: v(3),fc,pi2,ram,dum
	complex*16 :: cc,c11,c12,c21,c22,dumc
	real*8, allocatable :: gr(:,:,:),wr(:,:,:)
	complex*16, allocatable :: tmp(:),tmp1(:,:),
     &	                         sig_c(:,:,:),tmp2(:,:,:,:),
     &                           sig_tmp(:,:,:,:,:,:,:),
     &                           tmp3(:,:),g_c(:,:,:,:,:,:),
     &                           w_c(:,:,:),gc(:,:,:,:,:),
     &                           pcc(:,:,:,:,:),pc(:,:)
	pi2=pi+pi
	ram=1.d0/amega
	allocate(w_c(nplw_gw,maxval(ndim3_k),nr_full),stat=ierr)
      if(ierr/=0) call error_alloc('                       SIGC_II',
     &                             '                           W_C',
     &                 3,nplw_rm,maxval(ndim3_k),nr_full,idum,idum,
     &                 idum,idum,16)
	call w_q_g_r1(ind_tau,w_c)
	allocate(g_c(nbasmpw,maxval(ndim3_k),nr_full,nrel,2,nspin),
     &           stat=ierr)
      if(ierr/=0) call error_alloc('                       SIGC_II',
     &                             '                           G_C',
     &                   6,nbasmpw,maxval(ndim3_k),nr_full,nrel,
     &                   2,nspin,idum,16)
	do ispin=1,nspin
	  do it=1,2
	    call g_k_g_r1(it,ind_tau,g_c(1,1,1,1,1,ispin),ispin)
	  enddo
	enddo
	allocate(sig_tmp(nplw_rm,maxval(ndim3_k),nr_full,nrel,nrel,
     &         2,nspin),stat=ierr)
      if(ierr/=0) call error_alloc('                       SIGC_II',
     &                             '                       SIG_TMP',
     &                   7,nplw_rm,maxval(ndim3_k),nr_full,nrel,nrel,
     &                   2,nspin,16)
c -------------------------------------------------------
	allocate(sig_c(nqdiv,ndim3_nr_red(me3_k+1),max_r_star))
	allocate(pcc(nr_red,max_r_star,nqdiv,nrel,nrel))
	allocate(pc(nqdiv,nr_red))
	allocate(wr(nqdiv,ndim3_nr_red(me3_k+1),max_r_star))
	allocate(tmp3(nr_full,max_r_star))
	if(irel<=1) allocate(gr(nqdiv,ndim3_nr_red(me3_k+1),max_r_star))
	if(irel==2) allocate(gc(nqdiv,ndim3_nr_red(me3_k+1),max_r_star,
     &	                      nrel,nrel))
c --------------------------------------------------------------------
	sig_tmp=(0.d0,0.d0)
      do ista=1,nr_red
        nsta=num_r_star(ista)
        call w_rr_r_r1_star(w_c,wr,nsta,ista)
	  do ispin=1,nspin
	    do it=1,2
	      if(irel<=1) then
	        call g_rr_r_r1_star(ista,nsta,g_c(1,1,1,1,it,ispin),gr,
     &	                          dumc)
	      else if(irel==2) then
	        call g_rr_r_r1_star(ista,nsta,g_c(1,1,1,1,it,ispin),dum,
     &	                          gc)
            endif
            do jrl=1,nrel
              do irl=1,nrel
                do ind_i=1,ndim3_nr_red(me3_k+1)
	            i=n3_mpi_nr_red(me3_k+1)+ind_i
	            if(irel<=1) then
	              do ind_j=1,nsta
	                do ir=1,nqdiv
	                  sig_c(ir,ind_i,ind_j)=-gr(ir,ind_i,ind_j)
     &	                                      *wr(ir,ind_i,ind_j)
	                enddo
	                call fft3(ndiv(1),ndiv(2),ndiv(3),1,
     &	                        sig_c(1,ind_i,ind_j),-1)
	              enddo
	            else if(irel==2) then
	              do ind_j=1,nsta
	                do ir=1,nqdiv
	                  sig_c(ir,ind_i,ind_j)=
     &	                  -gc(ir,ind_i,ind_j,irl,jrl)
     &                    *wr(ir,ind_i,ind_j)
	                enddo
	                call fft3(ndiv(1),ndiv(2),ndiv(3),1,
     &	                        sig_c(1,ind_i,ind_j),-1)
	              enddo
	            endif
	          enddo   !! over ind_i
	          pcc(:,:,:,irl,jrl)=(0.d0,0.d0)
	          do k=1,nqdiv
	            call zone1_number(pnt(1,k),rb0,ndiv,ind)
                  do ind_i=1,ndim3_nr_red(me3_k+1)
	              i=n3_mpi_nr_red(me3_k+1)+ind_i
	              do ind_j=1,nsta !! over r' inside the star IND_S
	                pcc(i,ind_j,k,irl,jrl)=sig_c(ind,ind_i,ind_j)
                    enddo
	            enddo
	          enddo
	          if(nproc_k/=1) call dgop(pcc(1,1,1,irl,jrl),
     &	                                 2*nr_red*max_r_star*nqdiv,
     &                                   '  +',comm_pnt)
	        enddo   !! over irl
	      enddo   !! over jrl
c --------------------------------------------------------------------
	      do ind_k=1,ndim3_k(me3_k+1)
	        k=n3_mpi_k(me3_k+1)+ind_k
              do jrl=1,nrel
                do irl=1,nrel
c ------- We fill in the first index completely -------------
	            do i=1,nr_full
	              i0=ir_ref(i)
	              ig=r_group(i)
	              ind=k_sym_0(k,ig)  ! IND=A^(-1)*K  in NQDIV-list
                    c11=(0.0d0,0.0d0)
                    c12=(0.0d0,0.0d0)
                    c21=(0.0d0,0.0d0)
                    c22=(0.0d0,0.0d0)
                    if(irel<=1) then
	                c11=(1.d0,0.d0)
                    else if(irl==1.and.jrl==1) then
                      c11=conjg(uj(1,ig))*uj(1,ig)
                      c21=conjg(uj(3,ig))*uj(1,ig)
                      c12=conjg(uj(1,ig))*uj(3,ig)
                      c22=conjg(uj(3,ig))*uj(3,ig)
                    else if(irl==2.and.jrl==1) then
                      c11=conjg(uj(2,ig))*uj(1,ig)
                      c21=conjg(uj(4,ig))*uj(1,ig)
                      c12=conjg(uj(2,ig))*uj(3,ig)
                      c22=conjg(uj(4,ig))*uj(3,ig)
                    else if(irl==1.and.jrl==2) then
                      c11=conjg(uj(1,ig))*uj(2,ig)
                      c21=conjg(uj(3,ig))*uj(2,ig)
                      c12=conjg(uj(1,ig))*uj(4,ig)
                      c22=conjg(uj(3,ig))*uj(4,ig)
                    else if(irl==2.and.jrl==2) then
                      c11=conjg(uj(2,ig))*uj(2,ig)
                      c21=conjg(uj(4,ig))*uj(2,ig)
                      c12=conjg(uj(2,ig))*uj(4,ig)
                      c22=conjg(uj(4,ig))*uj(4,ig)
                    endif
	              do ind_j=1,nsta
	                j=ind_r_star(ind_j,ista)
	                iv=itshift_r(:,i,ig)-itshift_r(:,j,ig)
	                j0=ip_r(j,ig)
	                j0s=ind_r_star_back(j0)
                      cc=(0.0d0,0.0d0)
                      if(irel<=1) then
	                  cc=pcc(i0,j0s,ind,1,1)
                      else if(irel==2) then
                        cc=c11*pcc(i0,j0s,ind,1,1)
     &                    +c21*pcc(i0,j0s,ind,2,1)
     &                    +c12*pcc(i0,j0s,ind,1,2)
     &                    +c22*pcc(i0,j0s,ind,2,2)
                      endif
	                tmp3(i,ind_j)=cc*phase_k_rr(k,iv(1),iv(2),iv(3))
                    enddo
c ----- Multiplication with the phase factors ------------------
                    v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &    	 	                             +r_pnt(3,i)*r0b0(:,3)
	              fc=pi2*dot_product(pnt(:,k),v)
	              cc=exp(dcmplx(0.d0,-fc))
	              tmp3(i,1:nsta)=tmp3(i,1:nsta)*cc
	            enddo
	            call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nsta,tmp3,-1)
                  do ind_j=1,nsta
	              tmp3(:,ind_j)=tmp3(:,ind_j)/nr_full
	            enddo
	            do i=1,nplw_r(k)
	              i0=indpw_r(i,k)
	              call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
	              do ind_j=1,nsta
	                j=ind_r_star(ind_j,ista)
	                sig_tmp(i,ind_k,j,irl,jrl,it,ispin)=
     &   	                  tmp3(ind,ind_j)
	              enddo
	            enddo
	          enddo   !!! over irl
	        enddo  !! over jrl
	      enddo   !! over ind_k
	    enddo   !! over it
	  enddo   !! over ispin
	enddo    !! over ista
c --------------------------------------------------------------- 
	deallocate(sig_c,wr,w_c,g_c,pcc,pc,tmp3)
	if(irel<=1) deallocate(gr)
	if(irel==2) deallocate(gc)
	allocate(tmp1(nplw_rm,nplw_rm))
	allocate(tmp2(nbndf,nbndf,2,maxval(ndim3_k)))
	allocate(tmp3(nplw_rm,nbndf))
	allocate(tmp(nr_full))
	do ispin=1,nspin
	  tmp2=(0.d0,0.d0)
	  do it=1,2
          do jrl=1,nrel
            do irl=1,nrel
	        do ind_k=1,ndim3_k(me3_k+1)
	          k=n3_mpi_k(me3_k+1)+ind_k
	          do j=1,nplw_r(k)
                  do i=1,nr_full
                    v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     & 		                                 +r_pnt(3,i)*r0b0(:,3)
	              fc=pi2*dot_product(pnt(:,k),v)
	              cc=exp(dcmplx(0.d0,fc))
	              tmp(i)=sig_tmp(j,ind_k,i,irl,jrl,it,ispin)*cc
	            enddo
	            call fft3(nrdiv(1),nrdiv(2),nrdiv(3),1,tmp,1)
	            do i=1,nplw_r(k)
	              i0=indpw_r(i,k)
	              call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
	              tmp1(j,i)=tmp(ind)/nr_full
	            enddo
	          enddo
c -------- Contribution to Sigma_Corr ---------------------------
	          nbnd=n_bnd(k,ispin)
	          call zgemm('n','n',nplw_r(k),nbnd,nplw_r(k),
     &	                   dcmplx(ram,0.d0),tmp1,nplw_rm,
     &                     sa(1,1,ind_k,jrl,ispin),nplw_rm,(0.d0,0.d0),
     &                     tmp3,nplw_rm)
	          call zgemm('c','n',nbnd,nbnd,nplw_r(k),(1.d0,0.d0),
     &                     sa(1,1,ind_k,irl,ispin),nplw_rm,tmp3,nplw_rm,
     &                     (1.d0,0.d0),tmp2(1,1,it,ind_k),nbndf)
c --------------------------------------------------------------
	        enddo   !! over ind_k
            enddo   !! over irl
          enddo   !! over jrl
	  enddo  !! over it
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
          call ferm_pack_tau(tmp2(1,1,1,ind_k),
     &	                     sig_c_tau(1,1,1,ind_tau,ind_k,ispin),
     &	                     n_bnd(k,ispin),nbndf,nbndf)
	  enddo
	enddo   !! over ispin
	deallocate(tmp1,sig_tmp,tmp2,tmp3,tmp)
      end
