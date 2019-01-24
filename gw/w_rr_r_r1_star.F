      subroutine w_rr_r_r1_star(w_c,wr,nsta,ind_s)
c	Calculates W(R;r;r') where R-all; r=1,nr_red+par; r'=1,nstar(r1) -
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: nsta,ind_s
	complex*16, intent(in) :: w_c(nplw_gw,maxval(ndim3_k),
     &                              nr_full)
	real*8, intent(out) :: wr(nqdiv,ndim3_nr_red(me3_k+1),max_r_star)
	integer :: ind_k,k,i,i0,ind_j,j,kk,k0,ig,j0,j0s,ind,i1,iv(3),
     &	         ind_i,irr
	real*8 :: v(3),fc,pi2
	complex*16 :: cc
	complex*16, allocatable :: tmp(:,:),v_i(:,:,:)
	allocate(tmp(nr_full,nr_full))
	pi2=pi+pi
	allocate(v_i(nqdiv,nr_red,max_r_star))
	v_i=(0.d0,0.d0)
	wr=0.d0
	do ind_k=1,ndim3_k(me3_k+1)
	  k=n3_mpi_k(me3_k+1)+ind_k
c ------------------------------------------------------
	  tmp=(0.d0,0.d0)
	  do i=1,nplwgw(k) 
	    i0=indpw_gw(i,k)
	    call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
	    do ind_j=1,nsta
	      j=ind_r_star(ind_j,ind_s)
	      tmp(ind,ind_j)=w_c(i,ind_k,j)
	    enddo
	  enddo
	  call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nsta,tmp,1)
c ----- Multiplication with the phase factor ------------------
        do i=1,nr_full
          v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &	 		                   +r_pnt(3,i)*r0b0(:,3)
	    fc=pi2*dot_product(pnt(:,k),v)
	    cc=exp(dcmplx(0.d0,fc))
	    tmp(i,1:nsta)=tmp(i,1:nsta)*cc
	  enddo
c ----------------------------------------------------------------
	  do kk=1,nqdiv
	    k0=i_kref(kk)
	    if(k0/=k) cycle
	    ig=k_group(kk)
	    call zone1_number(pnt(:,kk),rb0,ndiv,ind)
	    do ind_j=1,nsta
	      j=ind_r_star(ind_j,ind_s)
	      j0=ip_r(j,ig)
	      j0s=ind_r_star_back(j0)
	      do i=1,nr_red
		      i1=r_pnt_ind(i)
	        i0=ip_r(i1,ig)
	        iv=itshift_r(:,i1,ig)-itshift_r(:,j,ig)
	        v_i(ind,i,ind_j)=tmp(i0,j0s)
     &                        *phase_k_rr(kk,iv(1),iv(2),iv(3))
            enddo
	    enddo
	  enddo  !! over kk
	enddo   !! over ind_k
	if(nproc_k/=1) call dgop(v_i,2*nqdiv*nr_red*max_r_star,'  +',
     &	                       comm_pnt)
      do ind_j=1,nsta
	  j=ind_r_star(ind_j,ind_s)
        do ind_i=1,ndim3_nr_red(me3_k+1)
	    i=n3_mpi_nr_red(me3_k+1)+ind_i
	    call fft3(ndiv(1),ndiv(2),ndiv(3),1,v_i(1,i,ind_j),1)
	    do irr=1,nqdiv
	      wr(irr,ind_i,ind_j)=v_i(irr,i,ind_j)/dfloat(nqdiv)
	    enddo
	  enddo
	enddo
	deallocate(v_i,tmp)
      end
