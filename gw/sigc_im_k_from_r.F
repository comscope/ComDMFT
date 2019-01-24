      subroutine sigc_im_k_from_r(ispin,ind_tau,sig_mi)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ispin,ind_tau
	complex*16, intent(in) :: sig_mi(nr_red,nfun,2,nrel,
     &	                               ndim3_kk(me3_k+1))
	integer :: ind,i,i0,j,ia,ib,ind_k,k,ig,jatom,jsort,nbnd,i_t,nbas,
     &	         irl,ind_ir,ir,k0,k1
	real*8 :: v(3),fc,pi2,sqa
	complex*16 :: cc
	complex*16, allocatable :: tmp(:,:),tmp1(:,:),tmp2(:,:),
     &	                         tmp3(:,:,:),tmp4(:),
     &	                         p_c(:,:),pcc(:,:,:,:,:,:)
	allocate(p_c(nqdiv,nfun))
	allocate(pcc(nfun,nr_red,2,nrel,k_star_max,maxval(ndim3_k)))
	do i_t=1,2   !! here is over TAU and BETA-TAU
        do irl=1,nrel
	    do i0=1,nr_red
	      p_c=(0.d0,0.d0)
	      do ind_ir=1,ndim3_kk(me3_k+1)
	        ir=n3_mpi_kk(me3_k+1)+ind_ir
	        do j=1,nfun
	          p_c(ir,j)=sig_mi(i0,j,i_t,irl,ind_ir)
	        enddo
	      enddo
	      if(nproc_k/=1) call dgop(p_c,2*nqdiv*nfun,'  +',comm_pnt)
	      call fft3(ndiv(1),ndiv(2),ndiv(3),nfun,p_c,-1)
	      do ind_k=1,ndim3_k(me3_k+1)
	        k=n3_mpi_k(me3_k+1)+ind_k
     	        do k0=1,k_star(k)
	          k1=k_list(k0,k)
	          call zone1_number(pnt(1,k1),rb0,ndiv,ind)
	          pcc(:,i0,i_t,irl,k0,ind_k)=p_c(ind,:)
	        enddo
	      enddo
	    enddo
	  enddo
	enddo  !! over i_t
	deallocate(p_c)
	allocate(tmp(nr_full,nfun),tmp4(nfun))
	allocate(tmp1(nfun,nplw_rm))
	allocate(tmp2(nfun,nbndf))
	allocate(tmp3(nbndf,nbndf,2))
	pi2=pi+pi
	sqa=1.d0/sqrt(amega)
	do ind_k=1,ndim3_k(me3_k+1)
	  k=n3_mpi_k(me3_k+1)+ind_k
	  nbas=nbask(k)/nrel
	  tmp3=(0.d0,0.d0)
	  do i_t=1,2   !! here is over TAU and BETA-TAU
          do irl=1,nrel
c ------- We fill in the first index completely -------------
	      do i=1,nr_full
	        i0=ir_ref(i)
	        ig=r_group(i)
	        ind=k_sym_0(k,ig)  ! IND=A^(-1)*K  in NQDIV-list
	        ind=ind_k_star_back(ind) ! position in a given K-star
              if(irel<=1) then
                tmp4=pcc(:,i0,i_t,1,ind,ind_k)
              else if(irl==1) then
                tmp4=uj(1,ig)*pcc(:,i0,i_t,1,ind,ind_k)
     &              +uj(3,ig)*pcc(:,i0,i_t,2,ind,ind_k)
              else if(irl==2) then
                tmp4=uj(2,ig)*pcc(:,i0,i_t,1,ind,ind_k)
     &              +uj(4,ig)*pcc(:,i0,i_t,2,ind,ind_k)
              endif
		      call sym_val_left(ig,tmp4,tmp(i,:),1)
	        do jatom=1,natom
	          jsort=is(jatom)
	          v=tshift(:,jatom,ig)-tshift_r(:,i,ig)
	          fc=pi2*dot_product(pnt(:,k),v)
	          cc=dcmplx(cos(fc),sin(fc))
	          ia=io_lem(jatom)
	          ib=ia+lfunm(jsort)-1
	          tmp(i,ia:ib)=tmp(i,ia:ib)*cc
	        enddo
	      enddo
c ----- Multiplication with the phase factors ------------------
            do i=1,nr_full
              v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &	  	  	                       +r_pnt(3,i)*r0b0(:,3)
	        fc=2.d0*pi*dot_product(pnt(:,k),v)
	        cc=exp(dcmplx(0.d0,fc))
	        tmp(i,1:nfun)=tmp(i,1:nfun)*cc
	      enddo
c -----------------------------------------------------------
	      call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nfun,tmp,1)
	      tmp=tmp/nr_full
	      tmp1=(0.d0,0.d0)
	      do i=1,nplw_r(k)
	        i0=indpw_r(i,k)
	        call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
	        tmp1(1:nfun,i)=tmp(ind,1:nfun)
	      enddo
c -------- Contribution to Sigma_Corr ---------------------------
	      nbnd=n_bnd(k,ispin)
	      call zgemm('n','n',nfun,nbnd,nplw_r(k),dcmplx(sqa,0.d0),
     &	               tmp1,nfun,sa(1,1,ind_k,irl,ispin),nplw_rm,
     &                 (0.d0,0.d0),tmp2,nfun)
	      call zgemm('c','n',nbnd,nbnd,nfun,(1.d0,0.d0),
     &	               z_bnd(1,1,ind_k,ispin),nfun,tmp2,nfun,
     &                 (1.d0,0.d0),tmp3(1,1,i_t),nbndf)
          enddo   !! over irl
c --- Account of S_mt_int + conjg(S_int_mt) --------------------------
	    do j=1,nbnd
	      do i=1,j
	        tmp3(i,j,i_t)=tmp3(i,j,i_t)+conjg(tmp3(j,i,i_t))
	        if(i/=j) tmp3(j,i,i_t)=conjg(tmp3(i,j,i_t))
	      enddo
	    enddo
	  enddo   !! over i_t
        call ferm_pack_tau(tmp3,sig_c_tau(1,1,1,ind_tau,ind_k,ispin),
     &	                   nbnd,nbndf,nbndf)
c --------------------------------------------------------------
	enddo    !! over ind_k
	deallocate(tmp,tmp1,tmp2,tmp3,tmp4,pcc)
      end
