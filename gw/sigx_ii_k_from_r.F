      subroutine sigx_ii_k_from_r(ispin)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: ispin
	integer :: ind,i,i0,j,ia,ib,ic,j0,ind_k,k,ig,ind_s,ind_j,ir,j0s,
     &	         iv(3),nsta,jrl,irl,ibas0,jbas0,nbas,ind_i,k0,k1,len,
     &           ind_ind_s
	real*8 :: v(3),fc,pi2,dum
	complex*16 :: cc,c11,c12,c21,c22,dumc
	real*8, allocatable :: gr(:,:,:),vr(:,:,:)
	complex*16, allocatable :: tmp(:,:),tmp1(:,:),s(:,:),
     &	                         sig_c(:,:,:,:,:),sig_tmp(:,:,:,:,:),
     &                           tmp3(:,:,:,:),g_c(:,:,:,:),
     &                           v_c(:,:,:),gc(:,:,:,:,:),
     &                           pcc(:,:,:,:,:,:),pc(:,:)
	pi2=pi+pi
	allocate(v_c(nplw_gw,maxval(ndim3_k),nr_full))
	call v_q_g_r1(v_c)
	allocate(g_c(nbasmpw,maxval(ndim3_k),nr_full,nrel))
	call gx_k_g_r1(g_c,ispin)
	allocate(sig_tmp(nplw_rm,maxval(ndim3_k),nr_full,nrel,nrel))
c -------------------------------------------------------	
	allocate(tmp(nr_full,nr_full))	
	allocate(tmp3(nr_full,max_r_star,nrel,nrel))
	allocate(sig_c(nqdiv,ndim3_nr_red(me3_k+1),max_r_star,nrel,nrel))
	allocate(pcc(nr_red,max_r_star,k_star_max,maxval(ndim3_k),
     &             nrel,nrel))	
	allocate(pc(nqdiv,nr_red))
	allocate(vr(nqdiv,ndim3_nr_red(me3_k+1),max_r_star))
	if(irel<=1) allocate(gr(nqdiv,ndim3_nr_red(me3_k+1),max_r_star))
	if(irel==2) allocate(gc(nqdiv,ndim3_nr_red(me3_k+1),max_r_star,
     &	                      nrel,nrel))
	sig_tmp=(0.d0,0.d0)
	do ind_ind_s=1,ndimt_nr_red(me3_tau+1)
	  ind_s=nt_mpi_nr_red(me3_tau+1)+ind_ind_s
        nsta=num_r_star(ind_s)
        call w_rr_r_r1_star(v_c,vr,nsta,ind_s)
	  if(irel<=1) call g_rr_r_r1_star(ind_s,nsta,g_c,gr,dumc)
	  if(irel==2) call g_rr_r_r1_star(ind_s,nsta,g_c,dum,gc)
	  sig_c=(0.d0,0.d0)
	  if(irel<=1) then
	    do ind_j=1,nsta
            do ind_i=1,ndim3_nr_red(me3_k+1)
	        do ir=1,nqdiv
	          sig_c(ir,ind_i,ind_j,1,1)=-gr(ir,ind_i,ind_j)
     &                                   *vr(ir,ind_i,ind_j)
	        enddo
	        call fft3(ndiv(1),ndiv(2),ndiv(3),1,
     &	                sig_c(1,ind_i,ind_j,1,1),-1)
	      enddo
	    enddo
	  else if(irel==2) then
          do jrl=1,2
            do irl=1,2
	        do ind_j=1,nsta
                do ind_i=1,ndim3_nr_red(me3_k+1)
	            do ir=1,nqdiv
	              sig_c(ir,ind_i,ind_j,irl,jrl)=
     &	              -gc(ir,ind_i,ind_j,irl,jrl)
     &                *vr(ir,ind_i,ind_j)
	            enddo
	            call fft3(ndiv(1),ndiv(2),ndiv(3),1,
     &	                    sig_c(1,ind_i,ind_j,irl,jrl),-1)
	          enddo
	        enddo
	      enddo
	    enddo
	  endif
        do jrl=1,nrel
          do irl=1,nrel
	      do ind_j=1,nsta !! over r' inside the current star IND_S
	        pc=(0.d0,0.d0)
              do ind_i=1,ndim3_nr_red(me3_k+1)  !! irreducible r
	          i=n3_mpi_nr_red(me3_k+1)+ind_i
	          pc(:,i)=sig_c(:,ind_i,ind_j,irl,jrl)
	        enddo
	        if(nproc_k/=1) call dgop(pc,2*nqdiv*nr_red,'  +',
     &	                               comm_pnt)
	        do ind_k=1,ndim3_k(me3_k+1)
	          k=n3_mpi_k(me3_k+1)+ind_k
     	          do k0=1,k_star(k)
	            k1=k_list(k0,k)
	            call zone1_number(pnt(1,k1),rb0,ndiv,ind)
	            pcc(:,ind_j,k0,ind_k,irl,jrl)=pc(ind,:)
	          enddo
	        enddo
	      enddo
	    enddo
	  enddo
c ---------------------------------------------------------------------
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
c ------- We fill in TMP for the given K completely -------------
	    do i=1,nr_full
	      i0=ir_ref(i)
	      ig=r_group(i)
	      ind=k_sym_0(k,ig)  ! IND=A^(-1)*K  in NQDIV-list
	      ind=ind_k_star_back(ind) ! position in a given K-star
            do jrl=1,nrel
              do irl=1,nrel
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
	            j=ind_r_star(ind_j,ind_s)
	            iv=itshift_r(:,i,ig)-itshift_r(:,j,ig)
	            j0=ip_r(j,ig)
	            j0s=ind_r_star_back(j0)
                  cc=(0.0d0,0.0d0)
                  if(irel<=1) then
	              cc=pcc(i0,j0s,ind,ind_k,1,1)
                  else if(irel==2) then
                    cc=c11*pcc(i0,j0s,ind,ind_k,1,1)
     &                +c21*pcc(i0,j0s,ind,ind_k,2,1)
     &                +c12*pcc(i0,j0s,ind,ind_k,1,2)
     &                +c22*pcc(i0,j0s,ind,ind_k,2,2)
                  endif
	            tmp3(i,ind_j,irl,jrl)=cc*phase_k_rr(k,iv(1),iv(2),
     &	                                              iv(3))
                enddo
              enddo   !! over irl
	      enddo    !! over jrl
	    enddo    !! over i
          do jrl=1,nrel
            do irl=1,nrel
c ----- Multiplication with the phase factors ------------------
              do i=1,nr_full
                v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &	 	                             +r_pnt(3,i)*r0b0(:,3)
	          fc=pi2*dot_product(pnt(:,k),v)
	          cc=exp(dcmplx(0.d0,-fc))
	          tmp(i,1:nsta)=tmp3(i,1:nsta,irl,jrl)*cc
	        enddo
	        call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nsta,tmp,-1)
              do ind_j=1,nsta
	          tmp(:,ind_j)=tmp(:,ind_j)/nr_full
	        enddo
	        do i=1,nplw_r(k)
	          i0=indpw_r(i,k)
	          call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
	          do ind_j=1,nsta
	            j=ind_r_star(ind_j,ind_s)
	            sig_tmp(i,ind_k,j,irl,jrl)=tmp(ind,ind_j)
	          enddo
	        enddo
	      enddo   !!! over irl
	    enddo  !! over jrl
	  enddo   !! over ind_k
	enddo    !! over ind_s
	deallocate(sig_c,vr,v_c,g_c,tmp3,pcc,pc)
	if(irel<=1) deallocate(gr)
	if(irel==2) deallocate(gc)
	len=2*nplw_rm*maxval(ndim3_k)*nr_full*nrel**2
	if(nproc_tau/=1) call dgop(sig_tmp,len,'  +',comm_pnk)
	allocate(tmp1(nplw_rm,nplw_rm))
	allocate(s(nbndtm,nplw_rm))
	do ind_k=1,ndim3_k(me3_k+1)
	  k=n3_mpi_k(me3_k+1)+ind_k
	  nbas=nbask(k)/nrel
c ------ We form the overlap matrix -----------------------------
	  do j=1,nplw_r(k)
	    j0=indpw_r(j,k)
	    do i=1,nbas
	      i0=indgb(i,k)
            ia=igbs(1,j0)-igbs(1,i0)
            ib=igbs(2,j0)-igbs(2,i0)
            ic=igbs(3,j0)-igbs(3,i0)
            ind=indplw(ia,ib,ic)
	      if(complex_ro) then
	        s(i,j)=dcmplx(sovr(ind),sovi(ind))
	      else
	        s(i,j)=dcmplx(sovr(ind),0.d0)
	      endif
	    enddo
	  enddo
c ------------------------------------------------------------	    
        do jrl=1,nrel
          jbas0=nbas*(jrl-1)
          do irl=1,nrel
            ibas0=nbas*(irl-1)
            do i=1,nr_full
              v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     & 	  	                           +r_pnt(3,i)*r0b0(:,3)
	        fc=pi2*dot_product(pnt(:,k),v)
	        cc=exp(dcmplx(0.d0,fc))
	        tmp(i,1:nplw_r(k))=sig_tmp(1:nplw_r(k),ind_k,i,irl,jrl)*cc
	      enddo
	      call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nplw_r(k),tmp,1)
	      do i=1,nplw_r(k)
	        i0=indpw_r(i,k)
	        call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
	        do j=1,nplw_r(k)
	          tmp1(j,i)=tmp(ind,j)/nr_full
	        enddo
	      enddo
	      call zgemm('n','n',nbas,nplw_r(k),nplw_r(k),
     &	               (1.d0,0.d0),s,nbndtm,tmp1,nplw_rm,(0.d0,0.d0),
     &                 tmp,nr_full)
	      call zgemm('n','c',nbas,nbas,nplw_r(k),(1.d0,0.d0),
     &	               tmp,nr_full,s,nbndtm,(0.d0,0.d0),
     &                 pgg(ibas0+1,jbas0+1,ind_k),nbasmpw)
          enddo   !! over irl
        enddo   !! over jrl
	enddo  !! over ind_k
	pgg=pgg*amega
	deallocate(tmp,tmp1,s,sig_tmp)
      end
