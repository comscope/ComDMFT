      subroutine ferm_rs_from_ks_ii(nbm,gk,ak,gr,gc,n_bn,it)
c	-----------  G_II from G_full -----------------------
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: nbm,n_bn(npnt),it
	real*8, intent(in) :: gk(nbm,nbm,maxval(ndim3_k))
      complex*16, intent(in) :: ak(nbasmpw,nbm,npnt)
	real*8, intent(out) :: gr(nqdiv,ndim3_nr_red(me3_k+1),nr_full,2)
	complex*16, intent(out) :: gc(nqdiv,ndim3_nr_red(me3_k+1),nr_full,
     &	                            2,nrel,nrel)
	integer :: k,nbnd,i,i0,ind,j,j0,ind_k,ig,i1,ind_j,kk,ind_s,irr,
     &	         iv(3),nsta,j0s,nsize,nbas,irl,jrl,ibas0,ind_i
	real*8 :: v(3),fc,pi2
	complex*16 :: cc,c11,c12,c21,c22
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),tmp0(:,:),
     &                           tmp2(:,:),g_c(:,:,:,:),g_ii(:,:,:),
     &                           tmp3(:,:),g_c1(:,:,:,:,:)
      nsize=max(nbasmpw,max_r_star)
      allocate(tmp0(nbm,nbasmpw),tmp2(nbm,nbm),
     &         tmp3(nr_full,max_r_star))
	allocate(g_c(nbasmpw,maxval(ndim3_k),nr_full,nrel))
      allocate(tmp1(nr_full,nbasmpw))
      allocate(tmp(nbasmpw,nbasmpw))
	pi2=pi+pi
c ---------------------------------------------------------------------
	do ind_k=1,ndim3_k(me3_k+1)
	  k=n3_mpi_k(me3_k+1)+ind_k
	  nbnd=n_bn(k)
		call unpack_hermit(tmp2,gk(1,1,ind_k),nbnd,nbm,nbm,(0.d0,0.d0),
     &		               (1.d0,0.d0))
	  call zgemm('n','c',nbnd,nbask(k),nbnd,(1.d0,0.d0),tmp2,
     &	           nbm,ak(1,1,k),nbasmpw,(0.d0,0.d0),tmp0,nbm)
	  call zgemm('n','n',nbask(k),nbask(k),nbnd,(1.d0,0.d0),
     &	           ak(1,1,k),nbasmpw,tmp0,nbm,(0.d0,0.d0),tmp,nbasmpw)
        nbas=nbask(k)/nrel
        do irl=1,nrel
          ibas0=(irl-1)*nbas
c ----- Fill in the first zone for second index -----------------------
          tmp1=(0.d0,0.d0)
	    do i=1,nbas
	      i0=indgb(i,k)
	      call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
	      tmp1(ind,1:nbask(k))=tmp(1:nbask(k),ibas0+i)
	    enddo
	    call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nbask(k),tmp1,-1)
c ----- Multiplication with the phase factor ------------------
          do i=1,nr_full
            v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &		 	                     +r_pnt(3,i)*r0b0(:,3)
	      fc=-pi2*dot_product(pnt(:,k),v)
	      cc=exp(dcmplx(0.d0,fc))/amega
	      g_c(1:nbask(k),ind_k,i,irl)=tmp1(i,1:nbask(k))*cc
	    enddo
	  enddo
	enddo   !! over ind_k
	deallocate(tmp0,tmp2,tmp,tmp1)
c ---- Now we have full G(G,r') for given irreducible K point -------
      do ind_s=1,nr_red !! we work with the stars of r'
        nsta=num_r_star(ind_s)
	  allocate(g_c1(nr_full,nsta,npnt,nrel,nrel))
	  g_c1=(0.d0,0.d0)
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
          nbas=nbask(k)/nrel
	    do jrl=1,nrel
            do irl=1,nrel
              ibas0=(irl-1)*nbas
c ----- Fill in the first zone for first index ----------------------
	        do i=1,nbas
	          i0=indgb(i,k)
	          call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
	          do ind_j=1,nsta
	            j=ind_r_star(ind_j,ind_s)
	            g_c1(ind,ind_j,k,irl,jrl)=g_c(ibas0+i,ind_k,j,jrl)
	          enddo
	        enddo
	        call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nsta,
     &	                g_c1(1,1,k,irl,jrl),1)
c ----- Multiplication with the phase factor ------------------
              do i=1,nr_full
                v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &		  	                         +r_pnt(3,i)*r0b0(:,3)
	          fc=pi2*dot_product(pnt(:,k),v)
	          cc=exp(dcmplx(0.d0,fc))
	          g_c1(i,1:nsta,k,irl,jrl)=g_c1(i,1:nsta,k,irl,jrl)*cc
	        enddo
	      enddo   !! over irl
	    enddo   !! over jrl
	  enddo   !! over ind_k
	  if(nproc_k/=1) then
	    call dgop(g_c1,2*nr_full*nsta*npnt*nrel**2,'  +',
     &	            comm_pnt)
	  endif
	  if(irel<=1) allocate(g_ii(nqdiv,ndim3_nr_red(me3_k+1),nsta))
	  do kk=1,nqdiv
	    k=i_kref(kk)
	    ig=k_group(kk)
	    call zone1_number(pnt(:,kk),rb0,ndiv,ind)
	    do jrl=1,nrel
            do irl=1,nrel
              if(irel<=1) then
                do j=1,nsta
                  do i=1,nr_full
                    tmp3(i,j)=g_c1(i,j,k,1,1)
                  enddo
                enddo
              else if(irel==2) then
                c11=(0.0d0,0.0d0)
                c12=(0.0d0,0.0d0)
                c21=(0.0d0,0.0d0)
                c22=(0.0d0,0.0d0)
	          if(irl==1.and.jrl==1) then
	            c11=conjg(uj(1,ig))*uj(1,ig)
	            c12=conjg(uj(1,ig))*uj(3,ig)
	            c21=conjg(uj(3,ig))*uj(1,ig)
	            c22=conjg(uj(3,ig))*uj(3,ig)
	          else if(irl==2.and.jrl==1) then
	            c11=conjg(uj(2,ig))*uj(1,ig)
	            c12=conjg(uj(2,ig))*uj(3,ig)
	            c21=conjg(uj(4,ig))*uj(1,ig)
	            c22=conjg(uj(4,ig))*uj(3,ig)
	          else if(irl==1.and.jrl==2) then
	            c11=conjg(uj(1,ig))*uj(2,ig)
	            c12=conjg(uj(1,ig))*uj(4,ig)
	            c21=conjg(uj(3,ig))*uj(2,ig)
	            c22=conjg(uj(3,ig))*uj(4,ig)
	          else if(irl==2.and.jrl==2) then
	            c11=conjg(uj(2,ig))*uj(2,ig)
	            c12=conjg(uj(2,ig))*uj(4,ig)
	            c21=conjg(uj(4,ig))*uj(2,ig)
	            c22=conjg(uj(4,ig))*uj(4,ig)
	          endif
                do j=1,nsta
                  do i=1,nr_full
                    tmp3(i,j)=c11*g_c1(i,j,k,1,1)+c12*g_c1(i,j,k,1,2)
     &                       +c21*g_c1(i,j,k,2,1)+c22*g_c1(i,j,k,2,2)
                  enddo
                enddo
              endif
	        do ind_j=1,nsta
	          j=ind_r_star(ind_j,ind_s)
	          j0=ip_r(j,ig)
	          j0s=ind_r_star_back(j0)
                do ind_i=1,ndim3_nr_red(me3_k+1)
	            i=n3_mpi_nr_red(me3_k+1)+ind_i
		          i1=r_pnt_ind(i)
	            i0=ip_r(i1,ig)
	            iv=itshift_r(:,i1,ig)-itshift_r(:,j,ig)
	            if(irel<=1) then
	              g_ii(ind,ind_i,ind_j)=tmp3(i0,j0s)
     &                          *phase_k_rr(kk,iv(1),iv(2),iv(3))
	            else if(irel==2) then
	              gc(ind,ind_i,j,it,irl,jrl)=
     &	              tmp3(i0,j0s)*phase_k_rr(kk,iv(1),iv(2),iv(3))
                  endif
	          enddo   !! over ind_i
	        enddo   !! over ind_j
	      enddo   !! over irl
	    enddo   !! over jrl
	  enddo   !! over kk
	  deallocate(g_c1)
	  if(irel<=1) then
          do ind_j=1,nsta
	      j=ind_r_star(ind_j,ind_s)
            do ind_i=1,ndim3_nr_red(me3_k+1)
	        i=n3_mpi_nr_red(me3_k+1)+ind_i
	        call fft3(ndiv(1),ndiv(2),ndiv(3),1,g_ii(1,ind_i,ind_j),1)
              do irr=1,nqdiv
	          gr(irr,ind_i,j,it)=g_ii(irr,ind_i,ind_j)/dfloat(nqdiv)
	        enddo
	      enddo
	    enddo
	  else if(irel==2) then
	    do jrl=1,2
	      do irl=1,2
              do ind_j=1,nsta
	          j=ind_r_star(ind_j,ind_s)
                do ind_i=1,ndim3_nr_red(me3_k+1)
	            i=n3_mpi_nr_red(me3_k+1)+ind_i
	            call fft3(ndiv(1),ndiv(2),ndiv(3),1,
     &	                    gc(1,ind_i,j,it,irl,jrl),1)
                  do irr=1,nqdiv
	              gc(irr,ind_i,j,it,irl,jrl)=
     &	               gc(irr,ind_i,j,it,irl,jrl)/dfloat(nqdiv)
	            enddo
	          enddo
	        enddo
	      enddo
	    enddo
	  endif
        if(irel<=1) deallocate(g_ii)
	enddo  !! over ind_s
	deallocate(tmp3)
	deallocate(g_c)
      end
