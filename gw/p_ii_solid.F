      subroutine p_ii_solid(ind_tau,isp1,isp2,p)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: ind_tau,isp1,isp2
      real*8, intent(inout) :: p(n_pbtot,n_pbtot,ndim3_tau,
     &                           maxval(ndim3_k))
	integer :: ind_i,nsta,ind_j,i,ir,jrl,irl,ind_k,k,i0,ig,ind,j,
     &	         iv(3),j0,j0s,ista,ierr,idum
	real*8 :: v(3),fc,pi2,dum
	complex*16 :: cc,dumc
	real*8, allocatable :: gr1(:,:,:),gr2(:,:,:)
	complex*16, allocatable :: tmp(:),tmp1(:,:),p_c(:,:,:),
     &	                         p_tmp(:,:,:),g_c(:,:,:,:,:),
     &                           gc1(:,:,:,:,:),gc2(:,:,:,:,:),s(:,:),
     &                           pcc(:,:,:),tmp3(:,:)
	pi2=pi+pi
	allocate(g_c(nbasmpw,maxval(ndim3_k),nr_full,nrel,2),
     &           stat=ierr)
      if(ierr/=0) call error_alloc('                    P_II_SOLID',
     &                             '                           G_C',
     &                    5,nbasmpw,maxval(ndim3_k),nr_full,nrel,
     &                    2,idum,idum,16)
	call g_k_g_r1(1,ind_tau,g_c,isp1)
	call g_k_g_r1(2,ind_tau,g_c,isp2)
	allocate(p_tmp(nplw_rm,maxval(ndim3_k),nr_full),stat=ierr)
      if(ierr/=0) call error_alloc('                    P_II_SOLID',
     &                             '                         P_TMP',
     &                 3,nplw_rm,maxval(ndim3_k),nr_full,idum,idum,
     &                 idum,idum,16)
	allocate(p_c(nqdiv,ndim3_nr_red(me3_k+1),max_r_star))
	allocate(pcc(nr_red,max_r_star,nqdiv))
	allocate(tmp3(nr_full,max_r_star))
	if(irel<=1) then
	  allocate(gr1(nqdiv,ndim3_nr_red(me3_k+1),max_r_star))
	  allocate(gr2(nqdiv,ndim3_nr_red(me3_k+1),max_r_star))
	else if(irel==2) then
	  allocate(gc1(nqdiv,ndim3_nr_red(me3_k+1),max_r_star,nrel,nrel))
	  allocate(gc2(nqdiv,ndim3_nr_red(me3_k+1),max_r_star,nrel,nrel))
	endif
	p_tmp=(0.d0,0.d0)
c --- P_C (R,r,r') where R - all; r - irred; r' - current star
      do ista=1,nr_red
        nsta=num_r_star(ista)   !! size of the star
        p_c=(0.d0,0.d0)
	  if(irel<=1) then
	    call g_rr_r_r1_star(ista,nsta,g_c(1,1,1,1,1),gr1,dumc)
	    call g_rr_r_r1_star(ista,nsta,g_c(1,1,1,1,2),gr2,dumc)
	  else if(irel==2) then
	    call g_rr_r_r1_star(ista,nsta,g_c(1,1,1,1,1),dum,gc1)
	    call g_rr_r_r1_star(ista,nsta,g_c(1,1,1,1,2),dum,gc2)
	  endif
        do ind_i=1,ndim3_nr_red(me3_k+1)
	    i=n3_mpi_nr_red(me3_k+1)+ind_i !! irreducible r
	    do ind_j=1,nsta   !! over r' inside the current star ISTA
            if(irel<=1) then
	        do ir=1,nqdiv
	          p_c(ir,ind_i,ind_j)=p_c(ir,ind_i,ind_j)
     &	                           -gr1(ir,ind_i,ind_j)
     &	                           *gr2(ir,ind_i,ind_j)
	        enddo
            else if(irel==2) then
              do jrl=1,2
                do irl=1,2
	            do ir=1,nqdiv
	              p_c(ir,ind_i,ind_j)=p_c(ir,ind_i,ind_j)
     &	                           -gc1(ir,ind_i,ind_j,irl,jrl)
     &	                           *conjg(gc2(ir,ind_i,ind_j,irl,jrl))
                  enddo
	          enddo
	        enddo   !! over jrl
	      endif
	      call fft3(ndiv(1),ndiv(2),ndiv(3),1,p_c(1,ind_i,ind_j),-1)
	    enddo   !! over ind_j
	  enddo   !! over ind_i
	  pcc=(0.d0,0.d0)
	  do k=1,nqdiv
	    call zone1_number(pnt(1,k),rb0,ndiv,ind)
          do ind_i=1,ndim3_nr_red(me3_k+1)
	      i=n3_mpi_nr_red(me3_k+1)+ind_i !! irreducible r
	      do ind_j=1,nsta   !! over r' inside the current star ISTA
	        pcc(i,ind_j,k)=p_c(ind,ind_i,ind_j)
	      enddo
	    enddo
	  enddo
	  if(nproc_k/=1) call dgop(pcc,2*nr_red*max_r_star*nqdiv,'  +',
     &	                         comm_pnt)
c --------------------------------------------------------------------
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
c ----------- Restore all r --------------------------------          
	    do i=1,nr_full
	      i0=ir_ref(i)
	      ig=r_group(i)
	      ind=k_sym_0(k,ig)  ! IND=A^(-1)*K  in NQDIV-list
	      do ind_j=1,nsta
	        j=ind_r_star(ind_j,ista)
	        iv=itshift_r(:,i,ig)-itshift_r(:,j,ig)
	        j0=ip_r(j,ig)
	        j0s=ind_r_star_back(j0)
	        tmp3(i,ind_j)=pcc(i0,j0s,ind)
     &	                   *phase_k_rr(k,iv(1),iv(2),iv(3))
            enddo
c ----- Multiplication with the phase factors ------------------
            v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &    	 	                     +r_pnt(3,i)*r0b0(:,3)
	      fc=pi2*dot_product(pnt(:,k),v)
	      cc=exp(dcmplx(0.d0,-fc))
	      tmp3(i,1:nsta)=tmp3(i,1:nsta)*cc
	    enddo
	    call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nsta,tmp3,-1)
          do ind_j=1,nsta
	      tmp3(:,ind_j)=tmp3(:,ind_j)/nr_full
	    enddo
c ------------------------------------------------------------
	    do i=1,nplw_r(k)
	      i0=indpw_r(i,k)
	      call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
	      do ind_j=1,nsta
	        j=ind_r_star(ind_j,ista)
	        p_tmp(i,ind_k,j)=tmp3(ind,ind_j)
	      enddo
	    enddo
	  enddo   !! over ind_k
	enddo    !! over ista
	deallocate(p_c,pcc,tmp3)
	if(irel<=1) deallocate(gr1,gr2)
	if(irel==2) deallocate(gc1,gc2)
	allocate(tmp3(nr_full,nplw_rm))
	allocate(tmp1(nplw_rm,nplw_rm))
	do ind_k=1,ndim3_k(me3_k+1)
	  k=n3_mpi_k(me3_k+1)+ind_k
	  allocate(tmp(nr_full))
	  do j=1,nplw_r(k)
          do i=1,nr_full
            v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     & 	                             +r_pnt(3,i)*r0b0(:,3)
	      fc=pi2*dot_product(pnt(:,k),v)
	      cc=exp(dcmplx(0.d0,fc))
	      tmp(i)=p_tmp(j,ind_k,i)*cc
	    enddo
	    call fft3(nrdiv(1),nrdiv(2),nrdiv(3),1,tmp,1)
	    do i=1,nplw_r(k)
	      i0=indpw_r(i,k)
	      call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
	      tmp1(j,i)=tmp(ind)/nr_full
	    enddo
	  enddo
        deallocate(tmp)
c ------ We form the overlap matrix -----------------------------
	  allocate(s(nplw_r(k),nplwgw(k)))
        call s_overlap_m_b(s,nplw_r(k),nplwgw(k),indpw_r(1,k),
     &                     indpw_gw(1,k))
	  call zgemm('n','n',nplw_r(k),nplwgw(k),nplw_r(k),
     &	           (1.d0,0.d0),tmp1,nplw_rm,s,nplw_r(k),(0.d0,0.d0),
     &             tmp3,nr_full)
	  call zgemm('C','n',nplwgw(k),nplwgw(k),nplw_r(k),
     &	           (1.d0,0.d0),s,nplw_r(k),tmp3,nr_full,(0.d0,0.d0),
     &             tmp1,nplw_rm)
        deallocate(s)
        do j=1,nplwgw(k)
	    do i=1,j
	      p(n_pbmt+i,n_pbmt+j,ind_tau,ind_k)=
     &	      p(n_pbmt+i,n_pbmt+j,ind_tau,ind_k)+dreal(tmp1(i,j))
	        if(i/=j) p(n_pbmt+j,n_pbmt+i,ind_tau,ind_k)=
     &	        p(n_pbmt+j,n_pbmt+i,ind_tau,ind_k)+dimag(tmp1(i,j))
	    enddo
	  enddo
	enddo  !! over ind_k
	deallocate(tmp1,tmp3,p_tmp)
      end
