      subroutine g_rr_r_r1_star(ista,nsta,g_c,gr,gc)
c	-----------  G_II from G_full -----------------------
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ista,nsta
	complex*16, intent(in) :: g_c(nbasmpw,maxval(ndim3_k),
     &                          nr_full,nrel)
	real*8, intent(out) :: gr(nqdiv,ndim3_nr_red(me3_k+1),max_r_star)
	complex*16, intent(out) :: gc(nqdiv,ndim3_nr_red(me3_k+1),
     &	                            max_r_star,nrel,nrel)
	integer :: ind_k,k,nbas,jrl,irl,ibas0,i,i0,j,ind_j,ind_i,ind,k0,
     &	         kk,ig,j0,j0s,i1,iv(3),irr
	real*8 :: v(3),fc,pi2
	complex*16 :: cc,c11,c12,c21,c22
      complex*16, allocatable :: g_ii(:,:,:,:,:),tmp3(:,:),g_c1(:,:,:,:)
      allocate(tmp3(nr_full,max_r_star))
	pi2=pi+pi
	allocate(g_ii(nqdiv,nr_red,nsta,nrel,nrel))
	allocate(g_c1(nr_full,nsta,nrel,nrel))
	g_ii=(0.d0,0.d0)
	do ind_k=1,ndim3_k(me3_k+1)
	  k=n3_mpi_k(me3_k+1)+ind_k
        nbas=nbask(k)/nrel
        g_c1=(0.d0,0.d0)
	  do jrl=1,nrel
          do irl=1,nrel
            ibas0=(irl-1)*nbas
c ----- Fill in the first zone for first index ----------------------
	      do i=1,nbas
	        i0=indgb(i,k)
	        call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
	        do ind_j=1,nsta
	          j=ind_r_star(ind_j,ista)
	          g_c1(ind,ind_j,irl,jrl)=g_c(ibas0+i,ind_k,j,jrl)
	        enddo
	      enddo
	      call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nsta,
     &                g_c1(1,1,irl,jrl),1)
c ----- Multiplication with the phase factor ------------------
            do i=1,nr_full
              v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &     	                           +r_pnt(3,i)*r0b0(:,3)
	        fc=pi2*dot_product(pnt(:,k),v)
	        cc=exp(dcmplx(0.d0,fc))
	        g_c1(i,1:nsta,irl,jrl)=g_c1(i,1:nsta,irl,jrl)*cc
	      enddo
	    enddo   !! over irl
	  enddo   !! over jrl
	  do k0=1,k_star(k)
	    kk=k_list(k0,k)
	    ig=k_group(kk)
	    call zone1_number(pnt(:,kk),rb0,ndiv,ind)
	    do jrl=1,nrel
            do irl=1,nrel
              if(irel<=1) then
                do j=1,nsta
                  do i=1,nr_full
                    tmp3(i,j)=g_c1(i,j,1,1)
                  enddo
                enddo
              else if(irel==2) then
                c11=(0.0d0,0.0d0)
                c12=(0.0d0,0.0d0)
                c21=(0.0d0,0.0d0)
                c22=(0.0d0,0.0d0)
	          if(irl==1.and.jrl==1) then
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
                do j=1,nsta
                  do i=1,nr_full
                    tmp3(i,j)=c11*g_c1(i,j,1,1)+c12*g_c1(i,j,1,2)
     &                       +c21*g_c1(i,j,2,1)+c22*g_c1(i,j,2,2)
                  enddo
                enddo
              endif
	        do ind_j=1,nsta
	          j=ind_r_star(ind_j,ista)
	          j0=ip_r(j,ig)
	          j0s=ind_r_star_back(j0)
                do i=1,nr_red
	            i1=r_pnt_ind(i)
	            i0=ip_r(i1,ig)
	            iv=itshift_r(:,i1,ig)-itshift_r(:,j,ig)
	            g_ii(ind,i,ind_j,irl,jrl)=tmp3(i0,j0s)
     &	                           *phase_k_rr(kk,iv(1),iv(2),iv(3))
	          enddo   !! over i
	        enddo   !! over ind_j
	      enddo   !! over irl
	    enddo   !! over jrl
	  enddo   !! over k0
	enddo   !! over ind_k
	deallocate(g_c1)
	deallocate(tmp3)
	if(nproc_k/=1) call dgop(g_ii,2*nqdiv*nsta*nr_red*nrel**2,'  +',
     &	                       comm_pnt)
	if(irel<=1) then
        do ind_j=1,nsta
	    j=ind_r_star(ind_j,ista)
          do ind_i=1,ndim3_nr_red(me3_k+1)
	      i=n3_mpi_nr_red(me3_k+1)+ind_i
	      call fft3(ndiv(1),ndiv(2),ndiv(3),1,g_ii(1,i,ind_j,1,1),1)
	      do irr=1,nqdiv
	        gr(irr,ind_i,ind_j)=g_ii(irr,i,ind_j,1,1)/dfloat(nqdiv)
	      enddo
	    enddo
	  enddo
	else if(irel==2) then
	  do jrl=1,2
	    do irl=1,2
            do ind_j=1,nsta
	        j=ind_r_star(ind_j,ista)
              do ind_i=1,ndim3_nr_red(me3_k+1)
	          i=n3_mpi_nr_red(me3_k+1)+ind_i
	          call fft3(ndiv(1),ndiv(2),ndiv(3),1,
     &                    g_ii(1,i,ind_j,irl,jrl),1)
	          do irr=1,nqdiv
	            gc(irr,ind_i,ind_j,irl,jrl)=
     &               g_ii(irr,i,ind_j,irl,jrl)/dfloat(nqdiv)
	          enddo
	        enddo
	      enddo
	    enddo
	  enddo
	endif
      deallocate(g_ii)
      end
